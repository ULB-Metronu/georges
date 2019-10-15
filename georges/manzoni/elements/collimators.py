"""
TODO
"""
from .magnets import Drift as _Drift


class Collimator(_Drift):
    PARAMETERS = {
        'APERTYPE': (None, 'Aperture type (CIRCULAR, ELIPTIC or RECTANGULAR)'),
        'APERTURE': ([], ''),
    }


class CircularCollimator(Collimator):
    PARAMETERS = {
        'APERTYPE': ('CIRCULAR', 'Aperture type'),
    }


class EllipticalCollimator(Collimator):
    PARAMETERS = {
        'APERTYPE': ('ELLIPTICAL', 'Aperture type'),
    }


class RectangularCollimator(Collimator):
    PARAMETERS = {
        'APERTYPE': ('RECTANGULAR', 'Aperture type'),
    }
