import sys


class Material:
    def __str__(self):
        return self.__class__.__name__.lower()

    def __eq__(self, y):
        return str(self) == str(y)


_materials = ['Air',
              'Aluminum',
              'Beryllium',
              'Gold',
              'Graphite',
              'Hydrogen',
              'Lead',
              'Lexan',
              'Mylar',
              'Oxygen',
              'Tin',
              'Titanium',
              'Water',
              'Vacuum'
]

for m in _materials:
    setattr(sys.modules[__name__], m, type(m, (Material,), {})())
