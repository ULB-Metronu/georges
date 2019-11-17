from __future__ import annotations
from typing import Optional
import logging
import os
try:
    import uproot as _uproot
except (ImportError, ImportWarning):
    logging.error("Uproot is required for this module to work.")
import pandas as _pd

__all__ = [
    'BDSimOutput',
    'ReBDSimOutput',
]


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

    class Directory:
        def __init__(self, directory):
            """

            Args:
                directory:
            """
            def _build(n, c):
                if c.__name__.endswith('Directory'):
                    return Output.Directory(directory=self._directory[n])
                else:
                    return self._directory[n]

            self._directory = directory
            for name, cls in self._directory.iterclasses():
                setattr(self, name.decode('utf-8').split(';')[0].replace('-', '_'), _build(name, cls))

    class Tree:
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

    class Branch:
        def __init__(self, tree: Output.Tree):
            self._tree: Output.Tree = tree
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

    class Header(Output.Tree):
        pass

    class Beam(Output.Tree):

        def read_df(self) -> _pd.DataFrame:
            """

            Returns:

            """

            beam_df = _pd.DataFrame()
            # Particle Name
            for branch, name in {'Beam.GMAD::BeamBase.beamParticleName': 'particleName'
                                 }.items():
                data = [_.decode('utf-8') for _ in self._tree.array(branch=[branch])]
            beam_df[name] = data

            # Single value
            for branch, name in {'Beam.GMAD::BeamBase.beamEnergy': 'E0',
                                 'Beam.GMAD::BeamBase.X0': 'X0',
                                 'Beam.GMAD::BeamBase.Y0': 'Y0',
                                 'Beam.GMAD::BeamBase.Z0': 'Z0',
                                 'Beam.GMAD::BeamBase.S0': 'S0',
                                 'Beam.GMAD::BeamBase.Xp0': 'Xp0',
                                 'Beam.GMAD::BeamBase.Yp0': 'Yp0',
                                 'Beam.GMAD::BeamBase.Zp0': 'Zp0',
                                 'Beam.GMAD::BeamBase.T0': 'T0',
                                 'Beam.GMAD::BeamBase.tilt': 'tilt',
                                 'Beam.GMAD::BeamBase.sigmaT': 'sigmaT',
                                 'Beam.GMAD::BeamBase.sigmaE': 'sigmaE',
                                 'Beam.GMAD::BeamBase.betx': 'betx',
                                 'Beam.GMAD::BeamBase.bety': 'bety',
                                 'Beam.GMAD::BeamBase.alfx': 'alfx',
                                 'Beam.GMAD::BeamBase.alfy': 'alfy',
                                 'Beam.GMAD::BeamBase.emitx': 'emitx',
                                 'Beam.GMAD::BeamBase.emity': 'emity',
                                 'Beam.GMAD::BeamBase.dispx': 'dispx',
                                 'Beam.GMAD::BeamBase.dispy': 'dispy',
                                 'Beam.GMAD::BeamBase.sigmaX': 'sigmaX',
                                 'Beam.GMAD::BeamBase.sigmaXp': 'sigmaXp',
                                 'Beam.GMAD::BeamBase.sigmaY': 'sigmaY',
                                 'Beam.GMAD::BeamBase.sigma11': 'sigma11',
                                 'Beam.GMAD::BeamBase.sigma12': 'sigma12',
                                 'Beam.GMAD::BeamBase.sigma13': 'sigma13',
                                 'Beam.GMAD::BeamBase.sigma14': 'sigma14',
                                 'Beam.GMAD::BeamBase.sigma15': 'sigma15',
                                 'Beam.GMAD::BeamBase.sigma16': 'sigma16',
                                 'Beam.GMAD::BeamBase.sigma22': 'sigma22',
                                 'Beam.GMAD::BeamBase.sigma23': 'sigma23',
                                 'Beam.GMAD::BeamBase.sigma24': 'sigma24',
                                 'Beam.GMAD::BeamBase.sigma25': 'sigma25',
                                 'Beam.GMAD::BeamBase.sigma26': 'sigma26',
                                 'Beam.GMAD::BeamBase.sigma33': 'sigma33',
                                 'Beam.GMAD::BeamBase.sigma34': 'sigma34',
                                 'Beam.GMAD::BeamBase.sigma35': 'sigma35',
                                 'Beam.GMAD::BeamBase.sigma36': 'sigma36',
                                 'Beam.GMAD::BeamBase.sigma44': 'sigma44',
                                 'Beam.GMAD::BeamBase.sigma45': 'sigma45',
                                 'Beam.GMAD::BeamBase.sigma46': 'sigma46',
                                 'Beam.GMAD::BeamBase.sigma55': 'sigma55',
                                 'Beam.GMAD::BeamBase.sigma56': 'sigma56',
                                 'Beam.GMAD::BeamBase.sigma66': 'sigma66',
                                 }.items():
                beam_df[name] = self._tree.array(branch=[branch])
            self._df = beam_df
            return self._df

    class Options(Output.Tree):
        pass

    class Model(Output.Tree):

        def read_df(self) -> _pd.DataFrame:
            """

            Returns:

            """
            model_geometry_df = _pd.DataFrame()

            # Names and strings
            for branch, name in {'Model.componentName': 'NAME',
                                 'Model.componentType': 'TYPE',
                                 'Model.material': 'MATERIAL',
                                 'Model.beamPipeType': 'APERTYPE',
                                 }.items():
                data = [_.decode('utf-8') for _ in self._tree.array(branch=[branch])[0]]
                model_geometry_df[name] = data

            # Scalar
            for branch, name in {'Model.length': 'L',
                                 'Model.staS': 'AT_ENTRY',
                                 'Model.midS': 'AT_CENTER',
                                 'Model.endS': 'AT_EXIT',
                                 'Model.tilt': 'TILT',
                                 'Model.k1': 'K1',
                                 'Model.k2': 'K2',
                                 'Model.k3': 'K3',
                                 'Model.k4': 'K4',
                                 'Model.k5': 'K5',
                                 'Model.k6': 'K6',
                                 'Model.k7': 'K7',
                                 'Model.k8': 'K8',
                                 'Model.k9': 'K9',
                                 'Model.k10': 'K10',
                                 'Model.k11': 'K11',
                                 'Model.k12': 'K12',
                                 'Model.k1s': 'K1S',
                                 'Model.k2s': 'K2S',
                                 'Model.k3s': 'K3S',
                                 'Model.k4s': 'K4S',
                                 'Model.k5s': 'K5S',
                                 'Model.k6s': 'K6S',
                                 'Model.k7s': 'K7S',
                                 'Model.k8s': 'K8S',
                                 'Model.k9s': 'K9S',
                                 'Model.k10s': 'K10S',
                                 'Model.k11s': 'K11S',
                                 'Model.k12s': 'K12S',
                                 'Model.bField': 'B',
                                 'Model.e1': 'E1',
                                 'Model.e2': 'E2',
                                 'Model.hgap': 'HGAP'

                                 }.items():
                model_geometry_df[name] = self._tree.array(branch=[branch])[0]

            # Aperture
            for branch, name in {'Model.beamPipeAper1': 'APERTURE1',
                                 'Model.beamPipeAper2': 'APERTURE2',
                                 'Model.beamPipeAper3': 'APERTURE3',
                                 'Model.beamPipeAper4': 'APERTURE4'}.items():
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
            return self._df

    class Run(Output.Tree):
        def __init__(self, tree):
            super().__init__(tree)

        class Summary(Output.Branch):
            pass

    class Event(Output.Tree):
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
            self.histos = BDSimOutput.Event.Histos(self)
            self.primary = BDSimOutput.Event.Primary(self)
            self.samplers = {
                b.decode('utf-8').rstrip('.'):
                    BDSimOutput.Event.Sampler(b, tree=self)
                for b in self._tree.keys() if b.decode('utf-8') not in (
                    'Summary.',
                    'Primary.',
                    'PrimaryGlobal.',
                    'Eloss.',
                    'ElossVacuum.',
                    'ElossTunnel.',
                    'ElossWorld.',
                    'ElossWorldExit.',
                    'PrimaryFirstHit.',
                    'PrimaryLastHit.',
                    'ApertureImpacts.',
                    'Trajectory.',
                    'Histos.',
                    'Primary.',
                )}

        def read_df(self):
            pass

        class ELoss(Output.Branch):
            pass

        class ELossVacuum(Output.Branch):
            pass

        class ELossTunnel(Output.Branch):
            pass

        class ELossWorld(Output.Branch):
            pass

        class ELossWorldExit(Output.Branch):
            pass

        class PrimaryFirstHit(Output.Branch):
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
                self._df.columns = self._df.columns.str.replace(self.__class__.__name__ + '.', '')
                return self._df

        class PrimaryLastHit(PrimaryFirstHit):
            pass

        class ApertureImpacts(Output.Branch):
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
                self._df.columns = self._df.columns.str.replace(self.__class__.__name__ + '.', '')
                return self._df

        class Histos(Output.Branch):
            def read_df(self) -> _pd.DataFrame:
                self._tree._tree

        class Primary(Output.Branch):
            def read_df(self) -> _pd.DataFrame:
                self._df = self._tree._tree.pandas.df(branches=tuple(map(lambda _: f"{self.__class__.__name__}.{_}", (
                    'n',
                    'energy',
                    'x',
                    'y',
                    'z',
                    'xp',
                    'yp',
                    'zp',
                    'T',
                    'weight',
                    'partID',
                    'S',
                    'r',
                    'rp',
                    'phi',
                    'phip',
                    'theta',
                    'charge',
                    'kineticEnergy',
                    'mass',
                    'rigidity',
                    'isIon',
                    'ionA',
                    'ionZ',
                ))))
                self._df.columns = self._df.columns.str.replace(self.__class__.__name__ + '.', '')
                return self._df

        class Sampler(Output.Branch):
            def __init__(self, name, *args, **kwargs):
                self._name = name.decode('utf-8')
                super().__init__(**kwargs)

            def read_df(self) -> _pd.DataFrame:
                self._df = self._tree._tree.pandas.df(branches=tuple(map(lambda _: self._name + _, (
                    'n',
                    'S',
                    'x',
                    'y',
                    'z',
                    'xp',
                    'yp',
                    'zp',
                    'energy',
                    'T',
                    'weight',
                    'partID',
                    'parentID',
                    'trackID',
                    'modelID',
                    'turnNumber',
                    # 'r',
                    # 'rp',
                    # 'phi',
                    # 'phip',
                    # 'theta',
                    # 'charge',
                    # 'kineticEnergy',
                    # 'mass',
                    # 'rigidity',
                    # 'isIon',
                    # 'ionA',
                    # 'ionZ',
                    # 'nElectrons',
                ))))
                self._df.columns = self._df.columns.str.replace(self._name, '')
                return self._df


class ReBDSimOutput(Output):

    def __init__(self, filename: str = 'output.root', path: str = '.'):
        """
        Note: we purposedly expose the 'ModelTree' tree as "model" to keep the same name as for a regular BDSimOutput
        object. The "Model" directory, which in almost all cases will be unused is exposed as "model_dir".
        Args:
            filename:
            path:
        """
        super().__init__(filename, path)
        self.header = BDSimOutput.Header(tree=self._file['Header'])
        self.model = BDSimOutput.Model(tree=self._file['ModelTree'])
        self.beam = ReBDSimOutput.Directory(directory=self._file['Beam'])
        self.event = ReBDSimOutput.Directory(directory=self._file['Event'])
        self.run = ReBDSimOutput.Directory(directory=self._file['Run'])
        self.options = ReBDSimOutput.Directory(directory=self._file['Options'])
        self.model_dir = ReBDSimOutput.Directory(directory=self._file['Model'])
