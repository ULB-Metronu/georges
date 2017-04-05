import os.path
import pandas as pd
import numpy as np
import numpy.linalg as npl
import georges.beam as beam
import georges.madx.madx as madx
from georges.sequence_geometry import compute_derived_data

DEFAULT_EXT = 'csv'


def beamline_is_defined(method):
    def with_check_defined(self, *args, **kwargs):
        if not hasattr(self, '_Beamline__beamline'):
            print("Beamline is not defined.")
            return
        else:
            return method(self, *args, **kwargs)
    return with_check_defined


class BeamlineException(Exception):
    """Exception raised for errors in the Beamline module."""

    def __init__(self, m):
        self.message = m


class Beamline:
    """A beamline or accelerator model.

    The internal representation is essentially a set of pandas DataFrames.
    """

    def __init__(self, *args, **kwargs):
        # We need the kwargs first to parse the args correctly
        self.__path = kwargs.get('path', '.')
        self.__prefix = kwargs.get('prefix', '')
        self.__elements = kwargs.get('elements', None)
        self.__length = kwargs.get('length', 0)
        self.__strengths = None
        self.__beam = kwargs.get('beam', None)
        self.__flag_ptc = kwargs.get('ptc', False)
        self.__flag_g4 = kwargs.get('g4', False)
        self.__madx_input = None
        self.__beamline = None
        self.__context = {}

        # A beam must be defined
        if self.__beam is None:
            print("No beam defined!")

        # Some type inference to get the elements right
        # Elements as a file name
        if self.__elements and isinstance(self.__elements, str):
            self.__build_elements_from_file(self.__elements)
        # Elements as a list to be converted onto a DataFrame
        elif self.__elements and isinstance(self.__elements, list):
            self.__elements = pd.DataFrame(self.__elements)

        # Now let's process the args
        for i, arg in enumerate(args):
            # Some type inference to get the sequence right
            # Sequence from file name
            if isinstance(arg, str) and i == 0:
                self.__name = arg.upper()
                self.__build_from_files([arg])
            # Sequence from a list of files
            if isinstance(arg, list) and len(arg) > 0 and i == 0:
                self.__name = '_'.join(arg).upper()
                self.__build_from_files(arg)
            # Sequence from a pandas.DataFrame
            if isinstance(arg, pd.DataFrame) and i == 0:
                self.__name = getattr(arg, 'name', 'BEAMLINE')
                self.__beamline = arg

        if self.__beamline is None:
            print("No beamline defined.")
            return

        # Check before hand if the survey will need to be converted
        survey = False
        if self.__beamline.get(['AT_ENTRY', 'AT_CENTER', 'AT_EXIT']) is None:
            if self.__beamline.get(['X', 'Y']) is not None:
                survey = True

        # Expand elements onto MAD-X native elements
        if self.__elements is not None and self.__beamline is not None:
            self.__expand_elements_data()

        # Compute derived data until a fixed point sequence is reached
        self.__expand_sequence_data()

        # If the sequence is given as a survey, convert to s-positions
        if survey:
            self.__convert_survey_to_sequence()
            # Re-expand
            self.__expand_sequence_data()

        # Compute the sequence length if needed
        if self.__length == 0:
            self.__length = self.__beamline['AT_EXIT'].values[-1]

        # Flag to distinguish MAD-X generated elements from beamline elements
        self.__beamline['PHYSICAL'] = True

    @property
    def name(self):
        """The sequence name."""
        return self.__name

    @property
    def length(self):
        """The sequence length."""
        return self.__length

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
    def last_madx_input(self):
        """Provides the last flat input sent to MAD-X."""
        return self.__madx_input

    @property
    @beamline_is_defined
    def line(self):
        """The beamline representation."""
        self.__beamline.name = self.name
        self.__beamline.length = self.length
        return self.__beamline

    @property
    def context(self):
        """The current state of the beamline."""
        if self.__context.get('PARTICLE') is None and self.__beam is not None:
            self.__context['PARTICLE'] = self.__beam.particle
        if self.__context.get('PC') is None and self.__beam is not None:
            self.__context['PC'] = self.__beam.pc
        if self.__context.get('BETAREL') is None and self.__beam is not None:
            self.__context['BETAREL'] = self.__beam.betarel
        return self.__context

    @context.setter
    def context(self, c):
        self.__context = c

    def set(self, k, v):
        """Set a single variable in the context. Allows method chaining."""
        self.__context[k] = v
        return self

    @property
    @beamline_is_defined
    def twiss(self, **kwargs):
        """Compute the Twiss parameters of the beamline."""
        # Override the argument to Beamline
        self.__flag_ptc = kwargs.get('ptc', self.__flag_ptc)
        m = madx.Madx(beamline=self, path=self.__path, madx='/usr/local/bin/madx-dev')
        m.beam()
        m.twiss(ptc=self.__flag_ptc, centre=True)
        print(m.input)
        errors = m.run(self.context).fatals
        self.__madx_input = m.input
        if len(errors) > 0:
            print(m.input)
            print(errors)
            raise BeamlineException("MAD-X ended with fatal error.")
        madx_twiss = madx.read_madx_twiss(os.path.join(self.__path, 'twiss.outx'))
        self.__beamline = madx_twiss.merge(self.__beamline,
                                           left_index=True,
                                           right_index=True,
                                           how='outer',
                                           suffixes=('_TWISS', '')
                                           ).sort_values(by='S')
        return self.__beamline

    @property
    @beamline_is_defined
    def track(self, **kwargs):
        """Compute the distribution of the beam as it propagates through the beamline."""
        if kwargs.get('ptc', False):
            self.__flag_ptc = True
        m = madx.Madx(beamline=self, path=self.__path, madx='/usr/local/bin/madx-dev')
        m.beam()
        m.track(self.__beam.distribution, ptc=self.__flag_ptc)
        errors = m.run(self.context).fatals
        self.__madx_input = m.input
        if len(errors) > 0:
            print(m.input)
            print(errors)
            raise BeamlineException("MAD-X ended with fatal error.")
        if self.__flag_ptc:
            madx_track = madx.read_ptc_tracking(os.path.join(self.__path, ''))
        else:
            madx_track = madx.read_madx_tracking(os.path.join(self.__path, 'tracking.outxone')).dropna()
            madx_track['PY'] = pd.to_numeric(madx_track['PY'])
        madx_track['S'] = round(madx_track['S'], 8)
        tmp = madx_track.query('TURN == 1').groupby('S').apply(lambda g: beam.Beam(g[['X', 'PX', 'Y', 'PY', 'PT']]))
        self.__beamline['AT_CENTER_TRUNCATED'] = round(self.__beamline['AT_CENTER'], 8)
        if 'BEAM' in self.__beamline:
            self.__beamline.drop('BEAM', inplace=True, axis=1)
        self.__beamline = self.__beamline.merge(pd.DataFrame(tmp, columns=['BEAM']),
                                                left_on='AT_CENTER_TRUNCATED',
                                                right_index=True,
                                                how='left').sort_values(by='AT_CENTER')
        self.__beamline.drop('AT_CENTER_TRUNCATED', axis=1, inplace=True)
        self.__beamline.sort_values(by='S', inplace=True)
        return self.__beamline


    def __build_from_files(self, names):
        files = [os.path.splitext(n)[0] + '.' + (os.path.splitext(n)[1] or DEFAULT_EXT) for n in names]
        try:
            sequences = [
                pd.read_csv(os.path.join(self.__path, self.__prefix, f), index_col = 'NAME') for f in files
            ]
        except OSError:
            print("One of the file has not been found.")
            return
        self.__beamline = pd.concat(sequences)

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
        try:
            self.__elements = pd.read_csv(os.path.join(self.__path, file), index_col='NAME')
        except OSError:
            print("ERROR: File apparently not found")
            self.__elements = pd.DataFrame()

    def __expand_elements_data(self):
        if self.__elements is None:
            return
        self.__beamline = self.__beamline.merge(self.__elements,
                                                left_on='TYPE',
                                                right_index=True,
                                                how='left',
                                                suffixes=('', '_ELEMENT')
                                                )
        # Angle conversion
        self.__beamline['ANGLE'] = self.__beamline['ANGLE'] / 180.0 * np.pi

    def __convert_survey_to_sequence(self):
        s = self.__beamline
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
