# Define constants for indexing in numpy arrays
INDEX = {k: i for i, k in enumerate(
    (
        'CLASS_CODE',
        'LENGTH',
        'ANGLE',
        'E1',
        'E2',
        'K1',
        'K2',
        'K3',
        'K4',
        'APERTYPE_CODE',
        'APERTURE',
        'APERTURE_2',
        'HGAP',
        'FINT',
        'FE_A0',
        'FE_A1',
        'FE_A2',
        'FE_DPP',
        'FE_LOSS',
        'BRHO',
    )
)}

# For convenience
INDEX_CLASS_CODE = INDEX['CLASS_CODE']
INDEX_LENGTH = INDEX['LENGTH']
INDEX_ANGLE = INDEX['ANGLE']
INDEX_E1 = INDEX['E1']
INDEX_E2 = INDEX['E2']
INDEX_K1 = INDEX['K1']
INDEX_K2 = INDEX['K2']
INDEX_K3 = INDEX['K3']
INDEX_K4 = INDEX['K4']
INDEX_APERTYPE_CODE = INDEX['APERTYPE_CODE']
INDEX_APERTURE = INDEX['APERTURE']
INDEX_APERTURE_2 = INDEX['APERTURE_2']
INDEX_HGAP = INDEX['HGAP']
INDEX_FINT = INDEX['FINT']
INDEX_FE_A0 = INDEX['FE_A0']
INDEX_FE_A1 = INDEX['FE_A1']
INDEX_FE_A2 = INDEX['FE_A2']
INDEX_FE_DPP = INDEX['FE_DPP']
INDEX_FE_LOSS = INDEX['FE_LOSS']
INDEX_BRHO = INDEX['BRHO']

# Define constants for elements types
CLASS_CODES = {k: i for i, k in enumerate(
    (
        'NONE',
        'DRIFT',
        'COLLIMATOR',
        'SBEND',
        'RBEND',
        'QUADRUPOLE',
        'SEXTUPOLE',
        'OCTUPOLE',
        'DECAPOLE',
        'MULTIPOLE',
        'ROTATION',
        'HKICKER',
        'VKICKER',
        'DEGRADER',
        'SCATTERER',
    )
)}

# Group elements in categories
CLASS_CODE_MATRIX = [
    CLASS_CODES['DRIFT'],
    CLASS_CODES['COLLIMATOR'],
    CLASS_CODES['QUADRUPOLE'],
    CLASS_CODES['SBEND'],
    CLASS_CODES['RBEND'],
    CLASS_CODES['ROTATION'],
]

CLASS_CODE_KICK = [
    CLASS_CODES['SEXTUPOLE'],
    CLASS_CODES['OCTUPOLE'],
    CLASS_CODES['DECAPOLE'],
    CLASS_CODES['MULTIPOLE'],
    CLASS_CODES['HKICKER'],
    CLASS_CODES['VKICKER'],
    CLASS_CODES['DEGRADER'],
    CLASS_CODES['SCATTERER'],

]

# Define constants for aperture types
APERTYPE_CODE_NONE = 0
APERTYPE_CODE_CIRCLE = 1
APERTYPE_CODE_RECTANGLE = 2

# Coordinates
X = 0
PX = 1
Y = 2
PY = 3
DPP = 4