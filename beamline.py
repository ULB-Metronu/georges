import os.path
import pandas as pd
import numpy as np
import numpy.linalg as npl
from .sequence_geometry import compute_derived_data


class BeamlineException(Exception):
    """Exception raised for errors in the Beamline module."""

    def __init__(self, m):
        self.message = m


class Beamline:
    """A beamline or accelerator model.

    The internal representation is essentially a pandas DataFrames.
    """

    def __init__(self, beamline, name=None, from_survey=False, with_expansion=True):
        """
        :param beamline: defines the beamline to be created. It can be
            - a single pandas Dataframe
            - a list that can be used to create a DataFrame
            - another beamline ('copy' operation)
        :param name: the name of the beamline to be created
        :param from_survey: indicates (True/False) is the beamline is to be converted from survey data
        """

        # Default values
        self.__length = 0
        self.__name = name
        self.__strengths = None
        self.__beamline = None
        self.__from_survey = from_survey

        # Process the beamline argument
        self.__create(beamline)

        # Verify that the beamline was created correctly
        if self.__beamline is None:
            raise BeamlineException("No beamline defined.")

        # Converts from survey data
        if self.__from_survey:
            self.__convert_angles_to_radians()
            self.__expand_sequence_data()
            self.__convert_survey_to_sequence()
            # Compute derived data until a fixed point sequence is reached,
            # again to reexpand

        # Compute derived data until a fixed point sequence is reached
        if with_expansion:
            self.__expand_sequence_data()

        # Compute the sequence length
        if self.__length == 0 and self.__beamline.get('AT_EXIT') is not None:
            self.__length = self.__beamline.get('AT_EXIT').max()

        # Beamline must be defined
        assert self.__length is not None
        assert self.__beamline is not None

    def __create(self, beamline):
        """Process the arguments of the initializer."""
        # Some type inference to get the sequence right
        # Sequence from a pandas.DataFrame
        if isinstance(beamline, pd.DataFrame):
            if 'NAME' in beamline.columns:
                self.__beamline = beamline.set_index('NAME') if beamline.index.names[0] is not 'NAME' else beamline
            else:
                self.__beamline = beamline
            if self.__name is None:
                self.__name = getattr(beamline, 'name', 'BEAMLINE')
            else:
                self.__beamline.name = self.__name
            if self.__beamline.size == 0:
                raise BeamlineException("Empty dataframe.")
        # Sequence from a list
        # Assume that a DataFrame can be created from the list
        if isinstance(beamline, list):
            self.__create(pd.DataFrame(beamline))
        # Sequence from another Beamline
        if isinstance(beamline, Beamline):
            self.__create(beamline.line)

    def __expand_sequence_data(self):
        """Apply sequence transformation until a fixed point is reached."""
        tmp = self.__beamline.apply(compute_derived_data, axis=1)
        tmp2 = tmp
        while True:
            tmp, tmp2 = tmp2, tmp.apply(compute_derived_data, axis=1)
            if tmp.equals(tmp2):
                break
        self.__beamline = tmp2

    def __convert_survey_to_sequence(self):
        s = self.__beamline
        if 'LENGTH' not in s:
            s['LENGTH'] = np.nan
        offset = s['ORBIT_LENGTH'][0] / 2.0
        if pd.isnull(offset):
            offset = 0
        self.__beamline['AT_CENTER'] = pd.DataFrame(
            npl.norm(
                [
                    s['X'].diff().fillna(0.0),
                    s['Y'].diff().fillna(0.0)
                ],
                axis=0
            ) - (
                s['LENGTH'].fillna(0.0) / 2.0 - s['ORBIT_LENGTH'].fillna(0.0) / 2.0
            ) + (
                s['LENGTH'].shift(1).fillna(0.0) / 2.0 - s['ORBIT_LENGTH'].shift(1).fillna(0.0) / 2.0
            )).cumsum() / 1000.0 + offset
        self.__converted_from_survey = True

    def __convert_angles_to_radians(self):
        # Angle conversion
        if 'ANGLE' in self.__beamline:
            self.__beamline['ANGLE'] *= np.pi / 180.0
        if 'ANGLE_ELEMENT' in self.__beamline:
            self.__beamline['ANGLE_ELEMENT'] *= np.pi / 180.0

    @property
    def name(self):
        """The sequence name."""
        return self.__name

    @property
    def length(self):
        """The sequence length."""
        return self.__length

    def add_extra_drift(self, extra):
        """Increase the sequence length by adding a drift"""
        self.__length += extra

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

    def add_drifts(self, with_pipe=True):
        line_with_drifts = pd.DataFrame()
        at_entry = 0

        def create_drift(name, length, at):
            s = pd.Series(
                {
                    'CLASS': 'DRIFT',
                    'TYPE': 'DRIFT',
                    'PIPE': with_pipe,
                    'LENGTH': length,
                    'AT_ENTRY': at,
                    'AT_CENTER': at + length / 2.0,
                    'AT_EXIT': at + length
                }
            )
            s.name = name
            return s

        for r in self.__beamline.iterrows():
            i = r[0]
            e = r[1]
            diff = e['AT_ENTRY'] - at_entry
            if diff == 0:
                line_with_drifts = line_with_drifts.append(e)
            else:
                line_with_drifts = line_with_drifts.append(create_drift(
                    name=f"DRIFT_{i}",
                    length=diff,
                    at=at_entry
                )).append(e)
                at_entry = e['AT_EXIT']

        return Beamline(line_with_drifts)
