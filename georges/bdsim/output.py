from __future__ import annotations
from typing import Optional, Any, List, Set
import multiprocessing
import concurrent.futures
from collections import UserDict
import logging
import glob
import os
try:
    import uproot as _uproot
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
    def __init__(self, filename: str = 'output.root', path: str = '.'):
        """
        Create a representation of a BDSIM output using uproot to read the root file.

        The root file is opened with uproot, so a valid path and filename must be provided.

        Args:
            filename: the name of the root file to read
            path: the path to the root file
        """
        self._executor = concurrent.futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count())
        self._files_pattern = os.path.join(path, filename)
        self._files: List[_uproot.rootio.ROOTDirectory] = []
        for file in glob.glob(self._files_pattern):
            f = self._executor.submit(_uproot.open, file)
            f.add_done_callback(lambda _: self._files.append(_.result()))

    def __getitem__(self, item):
        return self.files[0][item]

    @property
    def compression(self):
        """The compression algorithm used for the root file."""
        return self.files[0].compression

    @property
    def files(self) -> List[_uproot.rootio.ROOTDirectory]:
        """Return the expanded list of files."""
        if self._executor is not None:
            self._executor.shutdown(wait=True)
            self._executor = None
        return self._files

    class Directory:
        def __init__(self, output: Output, directory):
            """

            Args:
                directory:
            """
            def _build(n, c):
                if c.__name__.endswith('Directory'):
                    return Output.Directory(directory=self._directory[n])
                else:
                    return self._directory[n]

            self._output = output
            self._directory = directory
            for name, cls in self._directory.iterclasses():
                setattr(self, name.decode('utf-8').split(';')[0].replace('-', '_'), _build(name, cls))

        def __getitem__(self, item):
            return self._directory[item]

        @property
        def parent(self):
            return self._output

    class Tree:
        def __init__(self, output: Output, tree):
            """

            Args:
                tree: the associated uproot tree object
            """
            self._output = output
            self._df: Optional[_pd.DataFrame] = None
            self._np: Optional[_np.ndarray] = None
            self._executor = concurrent.futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count())
            self._tree = tree
            self._trees: List[Any] = []
            for file in output.files:
                f = self._executor.submit(file.get, tree)
                f.add_done_callback(lambda _: self._trees.append(_.result()))

        def __getitem__(self, item):
            if isinstance(item, list):
                return self.arrays(branches=item)
            else:
                return self.array(branch=item)

        def to_df(self) -> _pd.DataFrame:
            pass

        def to_np(self) -> _np.ndarray:
            pass

        def array(self, branch=None, **kwargs) -> _np.ndarray:
            """A proxy for the uproot method."""
            tmp = None
            for tree in self.trees:
                if tmp is None:
                    tmp = tree.array(branch=branch, **kwargs)
                else:
                    tmp = _np.concatenate([tmp, tree.array(branch=branch, **kwargs)])
            return tmp

        def arrays(self, branches=None, **kwargs):
            """A proxy for the uproot method."""
            return self._tree.arrays(branches=branches, **kwargs)

        @property
        def parent(self):
            return self._output

        @property
        def numentries(self):
            """Provides the number of entries in the tree (without reading the entire file)."""
            return self.trees[0].numentries

        @property
        def trees(self):
            """The associated uproot tree."""
            if self._executor is not None:
                self._executor.shutdown(wait=True)
                self._executor = None
            return self._trees

        @property
        def branches(self):
            return [b.decode('utf-8') for b in self.trees[0].keys()]

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

        def __init__(self, branch: str, tree: Output.Tree):
            self._branch: str = branch
            self._tree: Output.Tree = tree
            self._df: Optional[_pd.DataFrame] = None
            self._np: Optional[_np.ndarray] = None

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
                self._df = df
                #self._df.columns = self._df.columns.str.replace(self._branch, '')
            return self._df

        def to_np(self) -> _np.ndarray:
            pass

        def array(self, branch=None, **kwargs) -> _np.ndarray:
            """A proxy for the uproot method."""
            return self.parent.array(branch=self._branch + branch, **kwargs)

        def arrays(self, branches=None, **kwargs):
            """A proxy for the uproot method.
            TODO must be fixed
            """
            return self.parent.arrays(branches=[self._branch + b for b in branches], **kwargs)

        @property
        def parent(self) -> Output.Tree:
            return self._tree

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
            'beam',
            'options',
            'model',
            'run',
            'event',
        ):
            setattr(self,
                    item,
                    getattr(BDSimOutput, item.capitalize())(output=self, tree=item.capitalize())
                    )
            return getattr(self, item)

    class Header(Output.Tree):
        pass

    class Beam(Output.Tree):
        pass

    class Options(Output.Tree):
        pass

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
            for branch, name in {'Model.length': 'LENGTH',
                                 'Model.staS': 'AT_ENTRY',
                                 'Model.midS': 'AT_CENTER',
                                 'Model.endS': 'AT_EXIT',
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
                'T',
                'energy',
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
        # self.event = Output.Directory(directory=self.files[0]['Event'])
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
