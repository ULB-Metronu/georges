import os
import pandas as pd
from .beamline import Beamline

DEFAULT_EXTENSION = 'csv'


class BeamlineBuilderException(Exception):
    """Exception raised for errors in the BeamlineBuilder module."""

    def __init__(self, m):
        self.message = m


class BeamlineBuilder:
    def __init__(self, path='.', prefix='', elements=None):
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
        self.__beamline = pd.DataFrame()
        self.__name = ""
        self.__from_survey = False

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

    def add_from_survey_files(self, names, path=None, prefix=None):
        self.__from_survey = True
        return self.add_from_files(names, path, prefix)

    def add_from_file(self, file, path=None, prefix=None):
        return self.add_from_files([file], path, prefix)

    def add_from_survey_file(self, file):
        self.__from_survey = True
        return self.add_from_survey_files([file])

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
        self.__elements = pd.DataFrame(elements)
        return self

    def define_elements_from_file(self, file):
        file = os.path.splitext(file)[0] + '.' + (os.path.splitext(file)[1] or DEFAULT_EXTENSION)
        self.__elements = pd.read_csv(os.path.join(self.__path, file), index_col='NAME')
        return self

    def flatten(self, using='LENGTH', offset=0, drifts=0):
        if 'AT_CENTER' not in self.__beamline \
                and 'AT_EXIT' not in self.__beamline \
                and 'AT_ENTRY' not in self.__beamline:
            offset = offset
            at_entry_list = []
            for i, e in self.__beamline.iterrows():
                at_entry_list.append(offset)
                offset += e[using] + drifts
            self.__beamline['AT_ENTRY'] = at_entry_list
        return self

    def build(self, name=None, extra_length=0.0):
        if name is not None:
            self.__name = name
        if self.__elements is not None:
            self.__expand_elements_data()
        if 'PHYSICAL' not in self.__beamline:
            self.__beamline['PHYSICAL'] = True
        b = Beamline(self.__beamline, name=self.__name, from_survey=self.__from_survey)
        b.add_extra_drift(extra_length)
        return b

    def __expand_elements_data(self):
        self.__beamline = self.__beamline.merge(self.__elements,
                                                left_on='TYPE',
                                                right_index=True,
                                                how='left',
                                                suffixes=('', '_ELEMENT')
                                                )

    def add_element(self, e):
        e = pd.Series(e)
        e.name = e['NAME']
        self.__beamline = self.__beamline.append(e)
        return self

    # TODO
    def replicate(a, n=2):
        # return n*a
        pass

    def join(a, b):
        # return a + b
        pass


