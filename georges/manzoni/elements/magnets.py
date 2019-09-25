"""
TODO
"""
import numpy as _np
from ... import ureg as _ureg
from .elements import ManzoniElement as _ManzoniElement


class Marker(_ManzoniElement):
    """
    TODO
    """
    def propagate(self, beam: _np.ndarray, out: _np.ndarray = None):
        _np.copyto(out, beam, 'no')
        return out


class Drift(_ManzoniElement):
    PARAMETERS = {
        'L': (0.0 * _ureg.m, 'Drift length.'),
    }
    """Parameters of the element, with their default value and their descriptions."""


class Rotation(_ManzoniElement):
    PARAMETERS = {
        'ANGLE': (0.0 * _ureg.radian, 'Angle of rotation along the s-axis.')
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> list:
        return [
            self.ANGLE.m_as('radian')
        ]


class Quadrupole(_ManzoniElement):
    PARAMETERS = {
        'L': (0.0 * _ureg.m, 'Quadrupole length.'),
        'K1': (0.0 * _ureg.m**-2, 'Normalized gradient.'),
        'K1L': (0.0 * _ureg.m**-1, 'Normalized integrated gradient.'),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> list:
        return [
            self.L.m_as('m'), self.K1.m_as('m**-2'), 0.9
        ]


class Bend(_ManzoniElement):
    PARAMETERS = {
        'ANGLE': (0.0 * _ureg.radian, 'Bending angle.'),
        'K1': (0.0 * _ureg.m**-2, 'Quadrupolar normalized gradient.'),
        'K2': (0.0 * _ureg.m**-3, 'Sextupolar normalized gradient.'),
        'L': (0.0 * _ureg.m, 'Magnet length.'),
        'E1': (0.0 * _ureg.radian, 'Entrance face angle.'),
        'E2': (0.0 * _ureg.radian, 'Exirt face angle.'),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> list:
        return [
            self.L.m_as('m'), self.ANGLE.m_as('radian'), self.K1.m_as('m**-2'), self.K2.m_as('m**-3'), 0.9
        ]


class SBend(Bend):
    pass


class RBend(Bend):
    pass


class Solenoid(_ManzoniElement):
    pass


class Multipole(_ManzoniElement):
    pass


class Sextupole(_ManzoniElement):
    pass


class Octupole(_ManzoniElement):
    pass


class Decapole(_ManzoniElement):
    pass


class Dodecapole(_ManzoniElement):
    pass
