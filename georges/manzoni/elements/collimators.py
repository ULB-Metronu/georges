"""
TODO
"""
from ... import ureg as _ureg
from .elements import ManzoniElement as _ManzoniElement
from ..apertures import circular_aperture_check, rectangular_aperture_check


class Collimator(_ManzoniElement):
    PARAMETERS = {
        'APERTYPE:': ('', 'Aperture type'),
    }


class CircularCollimator(Collimator):
    PARAMETERS = {
        'APERTYPE': ('', 'Aperture type'),
    }


class RectangularCollimator(Collimator):
    PARAMETERS = {
        'APERTYPE': ('', 'Aperture type'),
    }
