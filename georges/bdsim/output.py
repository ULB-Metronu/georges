from __future__ import annotations
from typing import Optional
import logging
import os
try:
    import uproot as _uproot
except (ImportError, ImportWarning):
    logging.error("Uproot is required for this module to work.")
import pandas as _pd


class OutputType(type):
    pass


class Output(metaclass=OutputType):
    def __init__(self, filename: str = 'output.root', path: str = '.'):
        """

        Args:
            filename:
            path:
        """
        self._file = _uproot.open(os.path.join(path, filename))

    class _Tree:
        def __init__(self, tree):
            """

            Args:
                tree:
            """
            self._tree = tree
            self._df: Optional[_pd.DataFrame] = None

        def read_df(self):
            pass

        @property
        def numentries(self):
            return self._tree.numentries

        @property
        def df(self):
            if self._df is None:
                return self.read_df()
            return self._df

    class _Branch:
        def __init__(self, tree: Output._Tree):
            self._tree: Output._Tree = tree
            self._df: Optional[_pd.DataFrame] = None

        def read_df(self):
            pass

        @property
        def df(self):
            if self._df is None:
                return self.read_df()
            return self._df


class BDSimOutput(Output):
    def __init__(self, filename: str = 'output.root', path: str = '.'):
        """

        Args:
            filename:
            path:
        """
        super().__init__(filename, path)
        self.header = BDSimOutput.Header(tree=self._file['Header'])
        self.beam = BDSimOutput.Beam(tree=self._file['Beam'])
        self.options = BDSimOutput.Options(tree=self._file['Options'])
        self.model = BDSimOutput.Model(tree=self._file['Model'])
        self.run = BDSimOutput.Run(tree=self._file['Run'])
        self.event = BDSimOutput.Event(tree=self._file['Event'])

    class Header(Output._Tree):
        pass

    class Beam(Output._Tree):
        pass

    class Options(Output._Tree):
        pass

    class Model(Output._Tree):

        def read_df(self) -> _pd.DataFrame:
            """

            Returns:

            """
            model_geometry_df = _pd.DataFrame()

            # Names and strings
            for branch, name in {'Model.componentName': 'NAME',
                                 'Model.componentType': 'TYPE',
                                 'Model.material': 'MATERIAL'}.items():
                data = [_.decode('utf-8') for _ in self._tree.array(branch=[branch])[0]]
                model_geometry_df[name] = data

            # Scalar
            for branch, name in {'Model.length': 'LENGTH',
                                 'Model.staS': 'AT_ENTRY',
                                 'Model.midS': 'AT_CENTER',
                                 'Model.endS': 'AT_EXIT'}.items():
                model_geometry_df[name] = self._tree.array(branch=[branch])[0]

            # Vectors
            geometry_branches = {'Model.staPos': 'ENTRY_',
                                 'Model.midPos': 'CENTER_',
                                 'Model.endPos': 'EXIT_'}
            data = self._tree.pandas.df(branches=geometry_branches.keys(), flatten=True)
            for branch, name in geometry_branches.items():
                data.rename({f"{branch}.fX": f"{name}X", f"{branch}.fY": f"{name}Y", f"{branch}.fZ": f"{name}Z"},
                            axis='columns', inplace=True)

            # Concatenate
            self._df = _pd.concat([model_geometry_df, data.loc[0]], axis='columns', sort=False).set_index('NAME')

            self._df.columns = self._df.columns.str.lstrip(self.__class__.__name__ + '.')
            return self._df

    class Run(Output._Tree):
        def __init__(self, tree):
            super().__init__(tree)

        class Summary(Output._Branch):
            pass

    class Event(Output._Tree):
        def __init__(self, tree):
            super().__init__(tree)
            self.eloss = BDSimOutput.Event.ELoss(self)
            self.eloss_vacuum = BDSimOutput.Event.ELossVacuum(self)
            self.eloss_tunnel = BDSimOutput.Event.ELossTunnel(self)
            self.eloss_world = BDSimOutput.Event.ELossWorld(self)
            self.eloss_world_exit = BDSimOutput.Event.ELossWorldExit(self)
            self.primary_first_hit = BDSimOutput.Event.PrimaryFirstHit(self)
            self.primary_last_hit = BDSimOutput.Event.PrimaryLastHit(self)
            self.aperture_impacts = BDSimOutput.Event.ApertureImpacts(self)

        def read_df(self):
            pass

        class ELoss(Output._Branch):
            pass

        class ELossVacuum(Output._Branch):
            pass

        class ELossTunnel(Output._Branch):
            pass

        class ELossWorld(Output._Branch):
            pass

        class ELossWorldExit(Output._Branch):
            pass

        class PrimaryFirstHit(Output._Branch):
            def read_df(self) -> _pd.DataFrame:
                self._df = self._tree._tree.pandas.df(branches=tuple(map(lambda _: f"{self.__class__.__name__}.{_}", (
                    'n',
                    'S',
                    'weight',
                    #'partID',
                    #'trackID',
                    #'parentID',
                    'modelID',
                    'turn',
                    'x',
                    'y',
                    'z',
                    'X',
                    'Y',
                    'Z',
                    'T',
                    #'stepLength',
                    #'preStepKineticEnergy',
                    'storeTurn',
                    'storeLinks',
                    'storeModelID',
                    'storeLocal',
                    'storeGlobal',
                    'storeTime',
                    'storeStepLength',
                    'storePreStepKineticEnergy',
                ))))
                self._df.columns = self._df.columns.str.lstrip(self.__class__.__name__ + '.')
                return self._df

        class PrimaryLastHit(PrimaryFirstHit):
            pass

        class ApertureImpacts(Output._Branch):
            def read_df(self) -> _pd.DataFrame:
                self._df = self._tree._tree.pandas.df(branches=tuple(map(lambda _: f"{self.__class__.__name__}.{_}", (
                    'n',
                    'energy',
                    'S',
                    'weight',
                    'isPrimary',
                    'firstPrimaryImpact',
                    'partID',
                    'turn',
                    'x',
                    'y',
                    'xp',
                    'yp',
                    'T',
                    'kineticEnergy',
                    'isIon',
                    'ionA',
                    'ionZ',
                    'trackID',
                    'parentID',
                    'modelID',
                ))))
                self._df.columns = self._df.columns.str.lstrip(self.__class__.__name__ + '.')
                return self._df


class ReBDSimOutput(Output):
    pass
