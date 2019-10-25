"""
TODO
"""
from typing import Tuple
import numpy as _np
from ... import ureg as _ureg
from .elements import ManzoniElement as _ManzoniElement


class Marker(_ManzoniElement):
    """
    TODO
    """
    def propagate(self, beam_in: _np.ndarray, beam_out: _np.ndarray = None, parameters: list = None):
        _np.copyto(dst=beam_out, src=beam_in, casting='no')
        return beam_in, beam_out


class Drift(_ManzoniElement):
    PARAMETERS = {
        'L': (0.0 * _ureg.m, 'Drift length.'),
    }
    """Parameters of the element, with their default value and their descriptions."""

    def propagate(self,
                  beam_in: _np.ndarray,
                  beam_out: _np.ndarray = None,
                  global_parameters: list = None,
                  ) -> Tuple[_np.ndarray, _np.ndarray]:
        if self.L.magnitude == 0:
            _np.copyto(dst=beam_out, src=beam_in, casting='no')
            return beam_in, beam_out
        else:
            return self.integrator.propagate(self, beam_in, beam_out, global_parameters)

    @property
    def parameters(self) -> list:
        return list(map(float, [
            self.L.m_as('meter'),
        ]))


class Rotation(_ManzoniElement):
    PARAMETERS = {
        'ANGLE': (0.0 * _ureg.radian, 'Angle of rotation along the s-axis.')
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> list:
        return list(map(float, [
            self.ANGLE.m_as('radian'),
        ]))


class Magnet(_ManzoniElement):
    PARAMETERS = {
        'APERTYPE': (None, 'Aperture type (CIRCULAR, ELIPTIC or RECTANGULAR)'),
        'APERTURE': ([], ''),
        'KINEMATICS': (None, 'Reference kinematics')
    }


class Quadrupole(Magnet):
    PARAMETERS = {
        'L': (0.0 * _ureg.m, 'Quadrupole length.'),
        'K1': (0.0 * _ureg.m**-2, 'Normalized gradient.'),
        'K1S': (0.0 * _ureg.m**-1, 'Normalized skew gradient.'),
        'TILT': (0.0 * _ureg.radian, 'Magnet tilt.'),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> list:
        return list(map(float, [
            self.L.m_as('m'),
            self.K1.m_as('m**-2'),
            self.K1S.m_as('m**-2'),
            self.TILT.m_as('radian'),
        ]))


class Bend(Magnet):
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
        return list(map(float, [
            self.L.m_as('m'),
            self.ANGLE.m_as('radian'),
            self.K1.m_as('m**-2'),
            self.K2.m_as('m**-3'),
            self.KINEMATICS.beta or 1.0
        ]))


class SBend(Bend):
    pass


class RBend(Bend):
    pass


class Solenoid(Magnet):
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
