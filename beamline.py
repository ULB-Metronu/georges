import os.path
import pandas as pd
import numpy as np
import numpy.linalg as npl
from .sequence_geometry import compute_derived_data

DEFAULT_EXT = 'csv'


class BeamlineException(Exception):
    """Exception raised for errors in the Beamline module."""

    def __init__(self, m):
        self.message = m


class Beamline:
    """A beamline or accelerator model.

    The internal representation is essentially a pandas DataFrames.
    """

    def __init__(self, *args, **kwargs):
        """
        :param args: defines the beamline to be created. It can be
            - a single pandas Dataframe containing an existing beamline
            - another beamline ('copy' operation)
            - a csv file or a list of csv files (looked up in path/prefix)
        :param kwargs: optional parameters include:
            - path: filepath to the root directory of the beamline description files (defaults to '.')
            - prefix: prefix for the beamline description files (defaults to '')
            - elements: elements description file (looked up in path/)

        """
        # We need the kwargs first to parse the args correctly
        self.__path = kwargs.get('path', '.')
        self.__prefix = kwargs.get('prefix', '')
        self.__elements = kwargs.get('elements', None)
        self.__survey = kwargs.get('survey', False)

        # Default values
        self.__length = 0
        self.__strengths = None
        self.__beamline = None
        self.__converted_from_survey = False

        # Process the beamline argument
        self.__process_args(args)
        if self.__beamline is None:
            raise BeamlineException("No beamline defined.")

        # Process the elements description
        if self.__elements is not None:
            self.__process_elements()
            self.__expand_elements_data()

        # Angle conversion
        if self.__survey:
            if 'ANGLE' in self.__beamline:
                self.__beamline['ANGLE'] *= np.pi / 180.0
            if 'ANGLE_ELEMENT' in self.__beamline:
                self.__beamline['ANGLE_ELEMENT'] *= np.pi / 180.0

        # Compute derived data until a fixed point sequence is reached
        self.__expand_sequence_data()

        # If the sequence is given as a survey, convert to s-positions
        if self.__survey:
            self.__convert_survey_to_sequence()
            self.__expand_sequence_data()

        # Compute the sequence length
        if self.__length == 0 and self.__beamline.get('AT_EXIT') is not None:
            self.__length = self.__beamline.get('AT_EXIT').max()

        # Beamline must be defined
        assert self.__length is not None
        assert self.__beamline is not None

    def __str__(self):
        return str(self.__beamline)

    def __repr__(self):
        return self.__beamline.to_html()

    def __process_args(self, args):
        """Process the arguments of the initializer."""
        if len(args) == 0 or len(args) > 1:
            raise BeamlineException("Single argument expected.")
        arg = args[0]
        # Some type inference to get the sequence right
        # Sequence from file name
        if isinstance(arg, str):
            self.__name = arg.upper()
            self.__build_from_files([arg])
        # Sequence from a list of files
        if isinstance(arg, list) and len(arg) > 0:
            self.__name = '_'.join(arg).upper()
            self.__build_from_files(arg)
        # Sequence from a pandas.DataFrame
        if isinstance(arg, pd.DataFrame):
            self.__name = getattr(arg, 'name', 'BEAMLINE')
            self.__beamline = arg
            self.__beamline['PHYSICAL'] = True
            if self.__beamline.size == 0:
                raise BeamlineException("Empty dataframe.")
        # Sequence from another Beamline
        if isinstance(arg, Beamline):
            self.__name = arg.name
            self.__beamline = arg.line

    def __process_elements(self):
        """Process the elements description argument."""
        # Some type inference to get the elements right
        # Elements as a file name
        if isinstance(self.__elements, str):
            self.__build_elements_from_file(self.__elements)
        # Elements as a list to be converted onto a DataFrame
        elif isinstance(self.__elements, list) and len(self.__elements) > 0:
            self.__elements = pd.DataFrame(self.__elements)
        elif not isinstance(self.__elements, pd.DataFrame):
            raise BeamlineException("Invalid data type for 'elements'.")

    @property
    def name(self):
        """The sequence name."""
        return self.__name

    @name.setter
    def name(self, n):
        self.__name = n

    @property
    def length(self):
        """The sequence length."""
        return self.__length

    @property
    def converted_from_survey(self):
        """True if the sequence has been converted from survey data."""
        return self.__converted_from_survey

    @property
    def strengths(self):
        """Strengths of the various elements of the beamline."""
        return self.__strengths

    @strengths.setter
    def strengths(self, strengths):
        if not strengths.index == self.__strengths:
            raise BeamlineException("Trying to set the strengths for an invalid elements list.")
        self.__strengths = strengths

    @property
    def elements(self):
        """Elements composing the beamline."""
        return self.__elements

    @property
    def line(self):
        """The beamline representation."""
        self.__beamline.name = self.name
        self.__beamline.length = self.length
        return self.__beamline

    @line.setter
    def line(self, line):
        self.__beamline = line
        self.__length = line.get('AT_EXIT').max()

    def __build_from_files(self, names):
        """Build a sequence DataFrame from a list of .csv files."""
        files = [os.path.splitext(n)[0] + '.' + (os.path.splitext(n)[1] or DEFAULT_EXT) for n in names]
        sequences = [
            pd.read_csv(os.path.join(self.__path, self.__prefix, f), index_col='NAME') for f in files
        ]
        self.__beamline = pd.concat(sequences)
        self.__beamline['PHYSICAL'] = True

    def __expand_sequence_data(self):
        """Apply sequence transformation until a fixed point is reached."""
        tmp = self.__beamline.apply(compute_derived_data, axis=1)
        tmp2 = tmp
        while True:
            tmp, tmp2 = tmp2, tmp.apply(compute_derived_data, axis=1)
            if tmp.equals(tmp2):
                break
        self.__beamline = tmp2

    def __build_elements_from_file(self, file):
        file = os.path.splitext(file)[0] + '.' + (os.path.splitext(file)[1] or DEFAULT_EXT)
        self.__elements = pd.read_csv(os.path.join(self.__path, file), index_col='NAME')

    def __expand_elements_data(self):
        self.__beamline = self.__beamline.merge(self.__elements,
                                                left_on='TYPE',
                                                right_index=True,
                                                how='left',
                                                suffixes=('', '_ELEMENT')
                                                )

    def __convert_survey_to_sequence(self):
        s = self.__beamline
        if 'LENGTH' not in s:
            s['LENGTH'] = np.nan
        offset = s['ORBIT_LENGTH'][0] / 2.0
        if pd.isnull(offset):
            offset = 0
        self.__beamline['AT_CENTER'] = pd.DataFrame(npl.norm([
            s['X'].diff().fillna(0.0),
            s['Y'].diff().fillna(0.0)
        ], axis=0) - (
            s['LENGTH'].fillna(0.0) / 2.0 - s['ORBIT_LENGTH'].fillna(0.0) / 2.0
        ) + (
            s['LENGTH'].shift(1).fillna(0.0) / 2.0 - s['ORBIT_LENGTH'].shift(1).fillna(0.0) / 2.0
        )).cumsum() / 1000.0 + offset
        self.__converted_from_survey = True

    def add_markers(self):
        s = self.__beamline
        markers = []

        def create_marker(r):
            if r['CLASS'] != 'MARKER' and r['CLASS'] != 'INSTRUMENT':
                m = pd.Series({
                    'TYPE': 'MARKER',
                    'CLASS': 'MARKER',
                    'NAME': r.name + '_IN',
                    'AT_CENTER': r['AT_ENTRY'],
                    'PHYSICAL': False
                })
                markers.append(m)
                m = pd.Series({
                    'TYPE': 'MARKER',
                    'CLASS': 'MARKER',
                    'NAME': r.name + '_OUT',
                    'AT_CENTER': r['AT_EXIT'],
                    'PHYSICAL': False
                })
                markers.append(m)
            return r

        s.apply(create_marker, axis=1)
        return Beamline(pd.concat([s, pd.DataFrame(markers).set_index('NAME')]).sort_values(by='AT_CENTER'))

    def to_thin(self, element, value):
        bl = self.__beamline
        bl.set_value(element, 'LENGTH', 0.0)
        bl.set_value(element, 'ORBIT_LENGTH', 0.0)
        bl.set_value(element, 'AT_ENTRY', bl.loc[element]['AT_CENTER'])
        bl.set_value(element, 'AT_EXIT', bl.loc[element]['AT_CENTER'])
        bl.set_value(element, 'CLASS', 'QUADRUPOLE')
       # bl.set_value(element, 'PLUG', 'K1')
       # bl.set_value(element, 'VALUE', "{{0, {} }}".format(value))
        bl.set_value(element, 'APERTYPE', np.nan)
        bl.set_value(element, 'APERTURE', np.nan)
