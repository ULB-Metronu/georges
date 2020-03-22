"""
A Pythonic way to analyze and work with Beam Delivery SIMulation (BDSIM) ROOT output files.

Design goals:
 - No dependency on (py)ROOT(py) is needed. The module uses `uproot` instead.
 - Enables and favors exploration of the ROOT files. No prior knowledge of the content should be required
 to discover the content.
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Optional, Any, List, Set
import multiprocessing
import concurrent.futures
from collections import UserDict
import logging
import glob
import os
try:
    import uproot as _uproot
    if TYPE_CHECKING:
        import uproot.source.compressed
except (ImportError, ImportWarning):
    logging.error("Uproot is required for this module to work.")
import numpy as _np
import pandas as _pd

__all__ = [
    'BDSimOutput',
    'ReBDSimOutput',
    'ReBDSimOpticsOutput',
]


class OutputType(type):
    """A generic type for BDSIM output classes."""
    pass


class Output(metaclass=OutputType):
    def __init__(self, filename: str = 'output.root', path: str = '.', *, open_file: bool = True):
        """
        Create a representation of a BDSIM output using uproot to read the root file.

        The root file is opened with uproot, so a valid path and filename must be provided.

        Args:
            filename: the name of the root file to read
            path: the path to the root file
            open_file: attempts to open the master file
        """
        self._file = os.path.join(path, filename)
        if open_file:
            self._directory: _uproot.rootio.ROOTDirectory = _uproot.open(self._file)

    @classmethod
    def from_directory(cls, directory: _uproot.rootio.ROOTDirectory) -> Output:
        """Create an `Output` object directly attached to an existing ROOT directory.

        Args:
            directory: an existing `uproot` `ROOTDirectory`
        """
        o = cls(open_master=False)
        o._file = None
        o._directory = directory
        return o

    def __getitem__(self, item: str):
        """Read an object from the ROOT file or directory by name."""
        return self._directory[item]

    @property
    def compression(self) -> _uproot.source.compressed.Compression:
        """The compression algorithm used for the root file or directory."""
        return self._directory.compression

    @property
    def directory(self) -> _uproot.rootio.ROOTDirectory:
        """Return the master directory attached to this output."""
        return self._directory

    class Directory:
        def __init__(self, output: Output, directory: _uproot.rootio.ROOTDirectory):
            """
            A representation of a (nested) structure of ROOT directories.

            Args:
                output: the `Output` to which the directory structure is attached
                directory: the top-level ROOT directory
            """
            def _build(n, c):
                if c.__name__.endswith('Directory'):
                    return Output.Directory(output, directory=self._directory[n])
                else:
                    return self._directory[n]

            self._output: Output = output
            self._directory: _uproot.rootio.ROOTDirectory = directory
            for name, cls in self._directory.iterclasses():
                setattr(self, name.decode('utf-8').split(';')[0].replace('-', '_'), _build(name, cls))

        def __getitem__(self, item):
            return self._directory[item]

        @property
        def compression(self) -> _uproot.source.compressed.Compression:
            """The compression algorithm used for the directory."""
            return self._directory.compression

        @property
        def parent(self) -> Output:
            """The parent Output to which the directory structure is attached."""
            return self._output

    class Tree:
        def __init__(self, output: Output, tree: str):
            """
            A representation of a ROOT TTree structure.

            Args:
                output: the `Output` to which the tree is attached
                tree: the tree name
            """
            self._output = output
            self._df: Optional[_pd.DataFrame] = None
            self._np: Optional[_np.ndarray] = None
            self._tree: uproot.rootio.TTree = output[tree]

        def __getitem__(self, item):
            try:
                return self._tree[item]
            except KeyError:
                return self._tree[item + '.']

        def array(self, branch=None, **kwargs) -> _np.ndarray:
            """A proxy for the `uproot` `arrat` method."""
            return self.tree.array(branch=branch, **kwargs)

        def arrays(self, branches=None, **kwargs):
            """A proxy for the `uproot` `arrays` method."""
            return self.tree.arrays(branches=branches, **kwargs)

        def pandas(self, branches=None, **kwargs):
            """A proxy for the `uproot` `pandas` method."""
            return self._tree.pandas.df(branches=branches, **kwargs)

        def to_df(self) -> _pd.DataFrame:
            pass

        def to_np(self) -> _np.ndarray:
            pass

        @property
        def parent(self):
            """The parent Output to which the tree structure is attached."""
            return self._output

        @property
        def tree(self) -> _uproot.rootio.TTree:
            """The associated uproot tree."""
            return self._tree

        @property
        def branches(self) -> List[str]:
            return [b.decode('utf-8') for b in self.tree.keys()]

        @property
        def numentries(self) -> int:
            """Provides the number of entries in the tree (without reading the entire file)."""
            return self.tree.numentries

        @property
        def df(self) -> _pd.DataFrame:
            if self._df is None:
                return self.to_df()
            return self._df

        @property
        def np(self) -> _np.ndarray:
            if self._np is None:
                return self.to_np()
            return self._np

    class Branch:
        LEAVES: Set[str] = {}

        def __init__(self, tree: Output.Tree, branch: str):
            """
            A representation of a ROOT Branch.

            Args:
                tree: the `Tree` to which the branch is attached
                branch: the branch name
            """
            self._tree: Output.Tree = tree
            self._branch: str = tree[branch]
            self._df: Optional[_pd.DataFrame] = None
            self._np: Optional[_np.ndarray] = None

        def __getitem__(self, item):
            return self._branch[item]

        def array(self, branch=None, **kwargs) -> _np.ndarray:
            """A proxy for the `uproot` `array` method."""
            return self.parent.array(branch=self.branch.name.decode('utf-8') + branch, **kwargs)

        def arrays(self, branches=None, **kwargs):
            """A proxy for the uproot method.
            TODO must be fixed
            """
            return self.parent.arrays(branches=[self.branch.name + b for b in branches], **kwargs)

        def pandas(self, branches=None, **kwargs):
            """A proxy for the uproot method.
            TODO must be fixed
            """
            if branches is None:
                branches = self.leaves
            return self.parent.tree.pandas.df(branches=branches, **kwargs)

        def to_df(self) -> _pd.DataFrame:
            if self._df is None:
                df = _pd.DataFrame()
                for tree in self.parent.trees:
                    if len(self.LEAVES) == 0:
                        df = _pd.concat([
                            df,
                            tree.pandas.df()
                        ])
                    else:
                        df = _pd.concat([
                            df,
                            tree.pandas.df(
                                branches=tuple(map(lambda _: f"{self._branch}{_}", self.LEAVES))
                            )
                        ])
                df.columns = self.LEAVES
                self._df = df
            return self._df

        def to_np(self) -> _np.ndarray:
            pass

        @property
        def parent(self) -> Output.Tree:
            """The parent `Tree` to which the branch is attached."""
            return self._tree

        @property
        def branch(self) -> _uproot.rootio.TBranch:
            return self._branch

        @property
        def leaves(self) -> List:
            return self._branch.keys()

        @property
        def df(self) -> _pd.DataFrame:
            if self._df is None:
                return self.to_df()
            return self._df

        @property
        def np(self) -> _np.ndarray:
            if self._np is None:
                return self.to_np()
            return self._np


class BDSimOutput(Output):
    def __getattr__(self, item):
        if item in (
            'header',
            'geant4data',
            'beam',
            'options',
            'model',
            'run',
            'event',
        ):
            setattr(self,
                    item,
                    getattr(BDSimOutput, item.title())(output=self, tree=item.title())
                    )
            return getattr(self, item)

    class Header(Output.Tree):
        def __getattr__(self, b):
            if b in (
                    'header',
            ):
                setattr(self,
                        b,
                        getattr(BDSimOutput.Header, b.title())(branch=b.title() + '.', tree=self)
                        )
                return getattr(self, b)

        class Header(Output.Branch):
            LEAVES = {
                'bdsimVersion',
                'geant4Version',
                'rootVersion',
                'clhepVersion',
                'timeStamp',
                'fileType',
                'dataVersion',
                'doublePrecisionOutput',
                'analysedFiles',
                'combinedFiles',
                'nTrajectoryFilters',
                'trajectoryFilters',
            }

    class Beam(Output.Tree):
        def to_df(self) -> _pd.DataFrame:
            """

            Returns:

            """
            beam_df = _pd.Series()

            # Names and strings
            for branch, name in {'Beam.GMAD::BeamBase.particle': 'particleName',
                                 }.items():
                beam_df[name] = (self.trees[0].array(branch=[branch])[0]).decode('utf-8')

            # Single value
            for branch, name in {'Beam.GMAD::BeamBase.beamEnergy': 'E0',
                                 'Beam.GMAD::BeamBase.beamMomentum': 'P0',
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
                beam_df[name] = self._trees[0].array(branch=[branch])[0]

            self._df = _pd.DataFrame(beam_df).transpose()
            return self._df

    class Geant4Data(Output.Tree):
        ...

    class Options(Output.Tree):
        ...

    class Model(Output.Tree):
        def to_df(self) -> _pd.DataFrame:
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
                data = [_.decode('utf-8') for _ in self.trees[0].array(branch=[branch])[0]]
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
                                 'Model.eField': 'E',
                                 'Model.e1': 'E1',
                                 'Model.e2': 'E2',
                                 'Model.hgap': 'HGAP',
                                 'Model.fint': 'FINT',
                                 'Model.fintx': 'FINTX'
                                 }.items():
                model_geometry_df[name] = self.trees[0].array(branch=[branch])[0]

            # Aperture
            for branch, name in {'Model.beamPipeAper1': 'APERTURE1',
                                 'Model.beamPipeAper2': 'APERTURE2',
                                 'Model.beamPipeAper3': 'APERTURE3',
                                 'Model.beamPipeAper4': 'APERTURE4'}.items():
                model_geometry_df[name] = self.trees[0].array(branch=[branch])[0]

            # Vectors
            geometry_branches = {'Model.staPos': 'ENTRY_',
                                 'Model.midPos': 'CENTER_',
                                 'Model.endPos': 'EXIT_'}
            data = self.trees[0].pandas.df(branches=geometry_branches.keys(), flatten=True)
            for branch, name in geometry_branches.items():
                data.rename({f"{branch}.fX": f"{name}X", f"{branch}.fY": f"{name}Y", f"{branch}.fZ": f"{name}Z"},
                            axis='columns', inplace=True)

            # Concatenate
            self._df = _pd.concat([model_geometry_df, data.loc[0]], axis='columns', sort=False).set_index('NAME')

            return self._df

    class Run(Output.Tree):
        def __getattr__(self, item):
            if item in (
                    'summary',
            ):
                setattr(self,
                        item,
                        getattr(BDSimOutput.Run, item.capitalize())(branch='Summary.', tree=self)
                        )
                return getattr(self, item)

        class Summary(Output.Branch):
            pass

    class Event(Output.Tree):
        def __getattr__(self, item):
            if item in (
                    'eloss',
                    'eloss_vacuum',
                    'eloss_tunnel',
                    'eloss_world',
                    'eloss_world_exit',
                    'primary',
                    'primary_first_hit',
                    'primary_last_hit',
                    'aperture_impacts',
                    'histos'
            ):
                b = ''.join([i.capitalize() for i in item.split('_')])
                setattr(self,
                        item,
                        getattr(BDSimOutput.Event, b)(branch=b + '.', tree=self)
                        )
                return getattr(self, item)

            elif item == 'samplers':
                self.samplers = BDSimOutput.Event.Samplers({
                    b.decode('utf-8').rstrip('.'):
                        BDSimOutput.Event.Sampler(branch=b.decode('utf-8'), tree=self)
                    for b in self.trees[0].keys() if b.decode('utf-8') not in (
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
                    )})
                return self.samplers

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

        class Primary(Output.Branch):
            LEAVES = {
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
            }

        class PrimaryFirstHit(Output.Branch):
            LEAVES = {
                'n',
                'S',
                'weight',
                # 'partID',
                # 'trackID',
                # 'parentID',
                'modelID',
                'turn',
                'x',
                'y',
                'z',
                'X',
                'Y',
                'Z',
                'T',
                # 'stepLength',
                # 'preStepKineticEnergy',
                'storeTurn',
                'storeLinks',
                'storeModelID',
                'storeLocal',
                'storeGlobal',
                'storeTime',
                'storeStepLength',
                'storePreStepKineticEnergy',
            }

        class PrimaryLastHit(PrimaryFirstHit):
            pass

        class ApertureImpacts(Output.Branch):
            LEAVES = {
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
            }

        class Histos(Output.Branch):
            def read_df(self) -> _pd.DataFrame:
                pass

        class Samplers(UserDict):
            def compute_optics(self, samplers: Optional[List[str]] = None):
                return _pd.DataFrame(
                    [sampler.compute_optics() for sampler in self.data.values()]
                )

            def to_df(self, samplers: Optional[List[str]] = None, columns: Optional[List[str]] = None) -> _pd.DataFrame:
                pass

            def to_np(self, samplers: Optional[List[str]] = None, columns: Optional[List[str]] = None) -> _np.ndarray:
                pass

            @property
            def df(self) -> _pd.DataFrame:
                return self.to_df()

            @property
            def np(self) -> _np.ndarray:
                return self.to_np()

            @property
            def optics(self):
                if self._optics is None:
                    self._optics = self.compute_optics()
                return self._optics

        class Sampler(Output.Branch):
            LEAVES = {
                'x',
                'xp',
                'y',
                'yp',
                'z',
                'zp',
                'T',
                'p',
                'energy',
                'p',
                'turnNumber',
                'parentID',
                'partID',
                'trackID',
                'weight',
                'n',
                'S',
            }

            def to_np(self,
                      turn_number: int = -1,
                      primary_only: bool = True,
                      ) -> _np.ndarray:
                df: _pd.DataFrame = self.to_df()
                data = df[['x', 'xp', 'y', 'yp', 'T', 'energy', 'n', 'S']].values
                validity = df[['turnNumber', 'parentID']].values
                if turn_number == - 1 and primary_only is False:
                    return data
                elif turn_number == -1 and primary_only is True:
                    return data[validity[1, :] == 0]
                elif primary_only is False:
                    return data[validity[0, :] == turn_number]
                else:
                    return data[_np.logical_and(validity[:, 1] == 0, validity[:, 0] == turn_number), :]

            def compute_optics(self):
                """

                Returns:

                """
                data = self.to_np(turn_number=1, primary_only=True)
                cv = _np.cov(data)
                eps_x = _np.sqrt(cv[0, 0] * cv[1, 1] - cv[0, 1] * cv[1, 0])
                eps_y = _np.sqrt(cv[2, 2] * cv[3, 3] - cv[2, 3] * cv[3, 2])
                return {
                    'BETA11': (cv[0, 0] - cv[0, 5] ** 2) / eps_x,
                    'BETA22': (cv[2, 2] - cv[2, 5] ** 2) / eps_y,
                    'ALPHA11': -cv[1, 1] / eps_x,
                    'ALPHA22': -cv[3, 3] / eps_y,
                    'DISP1': cv[0, 5] / 0.001,
                    'DISP2': 0.0,
                    'DISP3': cv[2, 5] / 0.001,
                    'DISP4': 0.0,
                    'EPSX': eps_x,
                    'EPXY': eps_y,
                    'n': data[:, -2].sum(),
                    'S': data[0, -1],
                }


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
        # self.header = BDSimOutput.Header(output=self, tree=self.files[0]['Header'])
        # self.model = BDSimOutput.Model(output=self, tree=self.files[0]['ModelTree'])
        # self.beam = Output.Directory(directory=self.files[0]['Beam'])
        self.event = Output.Directory(self, directory=self.files[0]['Event'])
        # self.run = Output.Directory(directory=self.files[0]['Run'])
        # self.options = Output.Directory(directory=self.files[0]['Options'])
        # self.model_dir = Output.Directory(directory=self.files[0]['Model'])


class ReBDSimOpticsOutput(ReBDSimOutput):
    def __getattr__(self, item):
        if item in (
            'optics',
        ):
            setattr(self,
                    item,
                    getattr(ReBDSimOpticsOutput, item.capitalize())(output=self, tree=item.capitalize())
                    )
            return getattr(self, item)
        else:
            return getattr(super(), item)

    class Optics(Output.Tree):
        def __getattr__(self, item):
            if item in (
                    'optics',
            ):
                b = ''.join([i.capitalize() for i in item.split('_')])
                setattr(self,
                        item,
                        getattr(ReBDSimOpticsOutput.Optics, b)(branch='', tree=self)
                        )
                return getattr(self, item)

        class Optics(Output.Branch):
            pass
