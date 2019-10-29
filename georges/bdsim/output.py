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


class BDSimOutput(Output):
    def __init__(self, filename: str = 'output.root', path: str = '.'):
        """

        Args:
            filename:
            path:
        """
        super().__init__(filename, path)

    @property
    def model(self):
        """"""
        return BDSimOutput.Model(self._file['Model'])

    @property
    def event(self):
        """"""
        return BDSimOutput.Event(self._file['Event'])

    class _Branch:
        def __init__(self, branch):
            """

            Args:
                branch:
            """
            self._branch = branch

    class Model(_Branch):

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

    class Event(_Branch):
        pass


class ReBDSimOutput(Output):
    pass






