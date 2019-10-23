import os
import uproot as _uproot
import pandas as _pd


class OutputType(type):
    pass


class Output(metaclass=OutputType):
    def __init__(self, filename: str = 'output.root', path: str = '.'):
        self._file = _uproot.open(os.path.join(path, filename))


class BDSimOutput(Output):
    def __init__(self, filename: str = 'output.root', path: str = '.'):
        super().__init__(filename, path)

    @property
    def model(self):
        return BDSimOutput.Model(self._file['Model'])

    @property
    def event(self):
        return BDSimOutput.Event(self._file['Event'])

    class _Branch:
        def __init__(self, branch):
            self._branch = branch

    class Model(_Branch):

        def extract_geometry(self) -> _pd.DataFrame:
            model_geometry_df = _pd.DataFrame(
                columns=['Name', 'Type', 'Material', 'Length', 'Start_s', 'Mid_s', 'End_s', 'Start_x', 'Start_y',
                         'Start_z', 'Mid_x', 'Mid_y', 'Mid_z', 'End_x', 'End_y', 'End_z']
            )

            # Branch Names corresponding to string elements
            Branches_dico_type0 = {'Model.componentName': 'Name', 'Model.componentType': 'Type',
                                   'Model.material': 'Material'}

            # Branch Names corresponding to general float values
            Branches_dico_type1 = {'Model.length': 'Length', 'Model.staS': 'Start_s', 'Model.midS': 'Mid_s',
                                   'Model.endS': 'End_s'}

            # Branch Names corresponding to 3D position float values
            Branches_dico_type2 = {'Model.staPos': 'Start_', 'Model.midPos': 'Mid_', 'Model.endPos': 'End_'}

            # Filling of the BDSIM geometry model DataFrame for the different types of information

            for branch in Branches_dico_type0.keys():
                df = self._branch.pandas.df(branches=[branch])
                branch_list = df[branch][0]
                branch_list = [x.decode('utf-8') for x in branch_list]

                model_geometry_df[Branches_dico_type0[branch]] = branch_list

            for branch in Branches_dico_type1.keys():
                df = self._branch.pandas.df(branches=[branch])
                branch_list = df[branch][0].tolist()

                model_geometry_df[Branches_dico_type1[branch]] = branch_list

            for branch in Branches_dico_type2.keys():
                df = self._branch.pandas.df(branches=[branch])

                # X position
                branch_list = df[branch + '.fX'][0].tolist()
                model_geometry_df[Branches_dico_type2[branch] + 'x'] = branch_list

                # Y position
                branch_list = df[branch + '.fY'][0].tolist()
                model_geometry_df[Branches_dico_type2[branch] + 'y'] = branch_list

                # Z position
                branch_list = df[branch + '.fZ'][0].tolist()
                model_geometry_df[Branches_dico_type2[branch] + 'z'] = branch_list

    class Event(_Branch):
        pass


class ReBDSimOutput(Output):
    pass






