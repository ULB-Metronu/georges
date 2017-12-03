import os
import numpy as np
import numpy.linalg as npl
import pandas as pd
from .beamline import Beamline

DEFAULT_EXTENSION = 'csv'


class BeamlineBuilderException(Exception):
    """Exception raised for errors in the BeamlineBuilder module."""

    def __init__(self, m):
        self.message = m


class BeamlineBuilder:
    def __init__(self, path='.', prefix='', elements=None, survey=False):
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
        self.__internal = None
        self.__path = path
        self.__prefix = prefix
        self.__elements = None
        self.__beamline = None
        self.__name = ""

    def add_from_files(self, names, path=None, prefix=None):
        if path is not None:
            self.__path = path
        if prefix is not None:
            self.__prefix = prefix
        self.__name = '_'.join(names).upper()
        files = [os.path.splitext(n)[0] + '.' + (os.path.splitext(n)[1] or DEFAULT_EXTENSION) for n in names]
        sequences = [
            pd.read_csv(os.path.join(self.__path, self.__prefix, f), index_col='NAME') for f in files
        ]
        self.__beamline = pd.concat(sequences)
        self.__beamline['PHYSICAL'] = True
        return self

    def add_from_file(self, file, path=None, prefix=None):
        return self.add_from_files([file], path, prefix)

    def add_from_survey_file(self):
        return self.add_from_file()

    def define_elements(self, e):
        """Process the elements description argument."""
        # Some type inference to get the elements right
        # Elements as a file name
        if isinstance(e, str):
            return self.define_elements_from_file(e)
        # Elements as a list to be converted onto a DataFrame
        elif isinstance(e, list) and len(e) > 0:
            return self.define_elements_from_list(e)
        elif not isinstance(self.__elements, pd.DataFrame):
            raise BeamlineBuilderException("Invalid data type for 'elements'.")

    def define_elements_from_list(self, elements):
        self.__elements = pd.DataFrame(self.__elements)
        return self

    def define_elements_from_file(self, file):
        file = os.path.splitext(file)[0] + '.' + (os.path.splitext(file)[1] or DEFAULT_EXTENSION)
        self.__elements = pd.read_csv(os.path.join(self.__path, file), index_col='NAME')
        return self

    def convert(self, from_survey=True):
        if from_survey is True:
            self.__convert_survey_to_sequence()
        return self

    def build(self):
        self.__expand_elements_data()
        return Beamline(self.__beamline, name=self.__name)

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

    def __expand_elements_data(self):
        self.__beamline = self.__beamline.merge(self.__elements,
                                                left_on='TYPE',
                                                right_index=True,
                                                how='left',
                                                suffixes=('', '_ELEMENT')
                                                )

    # TODO
    def replicate(a, n=2):
        # return n*a
        pass

    def join(a, b):
        # return a + b
        pass

    def add_bend(self, **kwargs):
        bend = {
            'name': kwargs.get('name', 'bend'),
            'angle': kwargs.get('angle', 0),
            'K1': kwargs.get('K1', 0),
        }
        self._beamline.append(bend)
        return self
