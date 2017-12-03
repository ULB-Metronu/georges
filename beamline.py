import os.path
import pandas as pd
import numpy as np
from .sequence_geometry import compute_derived_data


class BeamlineException(Exception):
    """Exception raised for errors in the Beamline module."""

    def __init__(self, m):
        self.message = m


class Beamline:
    """A beamline or accelerator model.

    The internal representation is essentially a pandas DataFrames.
    """

    def __init__(self, beamline, name=None):
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

        # Default values
        self.__length = 0
        self.__name = name
        self.__strengths = None
        self.__beamline = None

        # Process the beamline argument
        self.__create(beamline)

        # Verify that the beamline was created correctly
        if self.__beamline is None:
            raise BeamlineException("No beamline defined.")

        # Compute derived data until a fixed point sequence is reached
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

    def __create(self, beamline):
        """Process the arguments of the initializer."""
        # Some type inference to get the sequence right
        # Sequence from a pandas.DataFrame
        if isinstance(beamline, pd.DataFrame):
            if self.__name is None:
                self.__name = getattr(beamline, 'name', 'BEAMLINE')
            self.__beamline = beamline
            if self.__beamline.size == 0:
                raise BeamlineException("Empty dataframe.")
        # Sequence from another Beamline
        if isinstance(beamline, Beamline):
            self.__name = beamline.name
            self.__beamline = beamline.line

    def __expand_sequence_data(self):
        """Apply sequence transformation until a fixed point is reached."""
        tmp = self.__beamline.apply(compute_derived_data, axis=1)
        tmp2 = tmp
        while True:
            tmp, tmp2 = tmp2, tmp.apply(compute_derived_data, axis=1)
            if tmp.equals(tmp2):
                break
        self.__beamline = tmp2

    @property
    def name(self):
        """The sequence name."""
        return self.__name

    @property
    def length(self):
        """The sequence length."""
        return self.__length

    @property
    def line(self):
        """The beamline representation."""
        self.__beamline.name = self.name
        self.__beamline.length = self.length
        return self.__beamline

    @property
    def strengths(self):
        """Strengths of the various elements of the beamline."""
        return self.__strengths

    @strengths.setter
    def strengths(self, strengths):
        if not strengths.index == self.__strengths:
            raise BeamlineException("Trying to set the strengths for an invalid elements list.")
        self.__strengths = strengths

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

    def to_thin(self, element):
        bl = self.__beamline
        bl.at[element, 'LENGTH'] = 0.0
        bl.at[element, 'ORBIT_LENGTH'] = 0.0
        bl.at[element, 'AT_ENTRY'] = bl.loc[element]['AT_CENTER']
        bl.at[element, 'AT_EXIT'] = bl.loc[element]['AT_CENTER']
        bl.at[element, 'CLASS'] = 'QUADRUPOLE'
        bl.at[element, 'APERTYPE'] = np.nan
        bl.at[element, 'APERTURE'] = np.nan

    def convert_angles_to_radians(self):
        # Angle conversion
        if 'ANGLE' in self.__beamline:
            self.__beamline['ANGLE'] *= np.pi / 180.0
        if 'ANGLE_ELEMENT' in self.__beamline:
            self.__beamline['ANGLE_ELEMENT'] *= np.pi / 180.0
