import sys
import os
import pkg_resources


class Material:
    def __str__(self):
        return self.__class__.__name__.lower()

    def __eq__(self, y):
        return str(self) == str(y)


_DATA_PATH = pkg_resources.resource_filename('georges', 'fermi/')
_pstar_data_files = [os.path.splitext(f)[0] for f in os.listdir(os.path.join(_DATA_PATH, 'pstar')) if
                     os.path.isfile(os.path.join(_DATA_PATH, 'pstar', f)) and os.path.splitext(f)[1] == '.txt']
_srim_data_files = [os.path.splitext(f)[0] for f in os.listdir(os.path.join(_DATA_PATH, 'srim')) if
                    os.path.isfile(os.path.join(_DATA_PATH, 'srim', f)) and os.path.splitext(f)[1] == '.txt']

_materials = list(map(str.title, (_pstar_data_files + _srim_data_files + ['Vacuum'])))

for m in _materials:
    setattr(sys.modules[__name__], m, type(m, (Material,), {})())
