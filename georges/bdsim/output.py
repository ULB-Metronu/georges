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

    class _Branch:
        pass


class BDSimOutput(Output):
    def __init__(self, filename: str = 'output.root', path: str = '.'):
        """

        Args:
            filename:
            path:
        """
        super().__init__(filename, path)

    @property
    def header(self):
        """"""
        return BDSimOutput.Header(self._file['Header'])

    @property
    def beam(self):
        """"""
        return BDSimOutput.Beam(self._file['Beam'])

    @property
    def options(self):
        """"""
        return BDSimOutput.Options(self._file['Options'])

    @property
    def model(self):
        """"""
        return BDSimOutput.Model(self._file['Model'])

    @property
    def run(self):
        """"""
        return BDSimOutput.Run(self._file['Run'])

    @property
    def event(self):
        """"""
        return BDSimOutput.Event(self._file['Event'])

    class Header(Output._Tree):
        pass

    class Beam(Output._Tree):
        pass

    class Options(Output._Tree):
        pass

    class Model(Output._Tree):

        def extract_geometry(self) -> _pd.DataFrame:
            """

            Returns:

            """
            model_geometry_df = _pd.DataFrame()

            # Names and strings
            for branch, name in {'Model.componentName': 'NAME',
                                 'Model.componentType': 'TYPE',
                                 'Model.material': 'MATERIAL'}.items():
                data = [_.decode('utf-8') for _ in self._branch.array(branch=[branch])[0]]
                model_geometry_df[name] = data

            # Scalar
            for branch, name in {'Model.length': 'LENGTH',
                                 'Model.staS': 'AT_ENTRY',
                                 'Model.midS': 'AT_CENTER',
                                 'Model.endS': 'AT_EXIT'}.items():
                model_geometry_df[name] = self._branch.array(branch=[branch])[0]

            # Vectors
            geometry_branches = {'Model.staPos': 'ENTRY_',
                                 'Model.midPos': 'CENTER_',
                                 'Model.endPos': 'EXIT_'}
            data = self._branch.pandas.df(branches=geometry_branches.keys(), flatten=True)
            for branch, name in geometry_branches.items():
                data.rename({f"{branch}.fX": f"{name}X", f"{branch}.fY": f"{name}Y", f"{branch}.fZ": f"{name}Z"},
                            axis='columns', inplace=True)

            # Concatenate
            return _pd.concat([model_geometry_df, data.loc[0]], axis='columns', sort=False).set_index('NAME')

    class Run(Output._Tree):
        def __init__(self):
            super().__init__()

        class Summary(Output._Branch):
            def __init__(self):
                pass

    class Event(Output._Tree):
        def __init__(self):
            super().__init__()
            self._eloss = BDSimOutput.Event.ELoss()
            self._eloss_vacuum = BDSimOutput.Event.ELossVacuum()
            self._eloss_tunnel = BDSimOutput.Event.ELossTunnel()
            self._eloss_world = BDSimOutput.Event.ELossWorld()
            self._eloss_world_exit = BDSimOutput.Event.ELossWorldExit()

        @property
        def eloss(self):
            return self._eloss

        class ELoss(Output._Branch):
            def __init__(self):
                pass

        class ELossVacuum(Output._Branch):
            def __init__(self):
                pass

        class ELossTunnel(Output._Branch):
            pass

        class ELossWorld(Output._Branch):
            pass

        class ELossWorldExit(Output._Branch):
            pass

        class PrimaryFirstHit(Output._Branch):
            pass

        class PrimaryLastHit(Output._Branch):
            pass

        class ApertureImpacts(Output._Branch):
            pass


class ReBDSimOutput(Output):
    pass
