import functools
import itertools
import collections
import os
import pandas as pd
from .beamline import Beamline

DEFAULT_EXTENSION = 'csv'


def flatten(iterable, ltypes=collections.Iterable):
    """
    Helper function to flatten nested/irregular lists.
    Borrowed from https://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists
    :param iterable: iterable object to be flattened
    :param ltypes: object type
    """
    remainder = iter(iterable)
    while True:
        first = next(remainder)
        if isinstance(first, ltypes) and not isinstance(first, (str, bytes)):
            remainder = itertools.chain(first, remainder)
        else:
            yield first


class BeamlineBuilderException(Exception):
    """Exception raised for errors in the BeamlineBuilder module."""

    def __init__(self, m):
        self.message = m


class BeamlineBuilder:
    def __init__(self, path='.', prefix='', elements=None):
        """
        BeamlineBuilder class used for the generation of georges.Beamline objects.
        :param path: path for file lookup (default to '.')
        :param prefix: default prefix for file lookup (default to '')
        :param elements: name of the elements file (default to None)
        """
        self.__path = path
        self.__prefix = prefix
        self.__elements = None
        self.__beamline = pd.DataFrame()
        self.__name = ""
        self.__from_survey = False
        if elements is not None:
            self.add_element(elements)

    @property
    def line(self):
        """
        :return: beamline object being built ('snapshot' during the creation phase).
        """
        return self.__beamline

    def add_from_files(self, names, path=None, prefix=None, sep=','):
        if path is not None:
            self.__path = path
        if prefix is not None:
            self.__prefix = prefix
        self.__name = '_'.join(names).upper()
        files = [
            os.path.splitext(n)[0] + (os.path.splitext(n)[1] or f".{DEFAULT_EXTENSION}") for n in names
        ]
        sequences = [
            pd.read_csv(os.path.join(self.__path, self.__prefix, f), index_col='NAME', sep=sep) for f in files
        ]
        if len(sequences) >= 2:
            sequences[1]['AT_CENTER'] += sequences[0].iloc[-1]['AT_CENTER']
            if sequences[0].index[-1] == sequences[1].index[0]:
                self.__beamline = pd.concat([sequences[0][:-1], sequences[1][1:]])
            else:
                self.__beamline = pd.concat(sequences)
            return self
        self.__beamline = pd.concat(sequences)
        self.__beamline['PHYSICAL'] = True
        return self

    def add_from_survey_files(self, names, path=None, prefix=None, sep=','):
        self.__from_survey = True
        return self.add_from_files(names, path, prefix, sep=sep)

    def add_from_file(self, file, path=None, prefix=None, sep=','):
        return self.add_from_files([file], path, prefix, sep=sep)

    def add_from_survey_file(self, file, path=None, prefix=None, sep=','):
        self.__from_survey = True
        return self.add_from_survey_files([file], path, prefix, sep=sep)

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

    def define_elements_from_file(self, file, sep=','):
        file = os.path.splitext(file)[0] + '.' + (os.path.splitext(file)[1] or DEFAULT_EXTENSION)
        self.__elements = pd.read_csv(os.path.join(self.__path, file), index_col='NAME', sep=sep)
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

    def add(self, e):
        if isinstance(e, (dict, pd.Series)):
            return self.add_element(e)
        if isinstance(e, (list,)):
            return self.add_sequence(e)

    def add_element(self, e):
        e = pd.Series(e)
        e.name = e['NAME']
        self.__beamline = self.__beamline.append(e)
        return self

    def add_sequence(self, s):
        self.__beamline = self.__beamline.append(pd.DataFrame(s))
        return self

    @staticmethod
    def build_sequence(s, start_drift=0.0, end_drift=0.0, inter_drift=0.0):
        def drift(l):
            return {
                'TYPE': 'DRIFT',
                'LENGTH': l,
                'ORBIT_LENGTH': l,
            }
        from operator import add
        sequence = list()
        if start_drift != 0.0:
            sequence.append(drift(start_drift))
        if inter_drift != 0.0:
            for e in functools.reduce(add, [(e, drift(inter_drift)) for e in s])[:-1]:
                sequence.append(e)
        else:
            for e in s:
                sequence.append(e)
        if end_drift != 0.0:
            sequence.append(drift(end_drift))
        return list(flatten(sequence, ltypes=list))

    def flatten(self, using='LENGTH', offset=0):
        elements_names = {}
        at_entry_list = []
        for i, e in self.__beamline.iterrows():
            if elements_names.get(e['NAME']) is None:
                elements_names[e['NAME']] = 0
            else:
                elements_names[e['NAME']] += 1
            self.__beamline.loc[i, 'NAME'] = f"{self.__beamline.iloc[i]['NAME']}_{elements_names[e['NAME']]}"
            at_entry_list.append(offset)
            offset += e[using]
        self.__beamline['AT_ENTRY'] = at_entry_list
        if self.__beamline.iloc[-1]['TYPE'] is 'DRIFT':
            self.__beamline = self.__beamline.append(
                {
                    'NAME': 'END_BUILDER_MARKER',
                    'TYPE': 'MARKER',
                    'AT_ENTRY': at_entry_list[-1] + self.__beamline.iloc[-1][using]
                }
                , ignore_index=True
            )
        self.__beamline = self.__beamline.query("TYPE != 'DRIFT'")
        return self
