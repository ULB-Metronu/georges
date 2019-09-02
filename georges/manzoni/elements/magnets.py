"""
TODO
"""
from typing import Callable
import numpy as _np
from ... import ureg as _ureg
from .elements import ManzoniElement as _ManzoniElement
from ..kernels import drift5, batched_vector_matrix, batched_vector_matrix_tensor
from ..kernels import FirstOrderIntegrator, SecondOrderIntegrator


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

    def build_kernel(self):
        self._kernel = drift5
        self._kernel_arguments = [_np.array([self.L.m_as('m')])]


class Rotation(_ManzoniElement):
    PARAMETERS = {
        'ANGLE': (0.0 * _ureg.radian, 'Angle of rotation along the s-axis.')
    }
    """Parameters of the element, with their default value and their descriptions."""

    def build_kernel(self):
        self._kernel = batched_vector_matrix
        angle = self.ANGLE.m_as('radian')
        self._kernel_arguments = [_np.array([
            [_np.cos(angle), 0, -_np.sin(angle), 0, 0],
            [0, _np.cos(angle), 0, -_np.sin(angle), 0],
            [_np.sin(angle), 0, _np.cos(angle), 0, 0],
            [0, _np.sin(angle), 0, _np.cos(angle), 0],
            [0, 0, 0, 0, 1],
        ])]


class Quadrupole(_ManzoniElement):
    PARAMETERS = {
        'L': (0.0 * _ureg.m, 'Quadrupole length.'),
        'K1': (0.0 * _ureg.m**-2, 'Normalized gradient.'),
        'K1L': (0.0 * _ureg.m**-1, 'Normalized integrated gradient.'),
    }
    """Parameters of the element, with their default value and their descriptions."""

    def build_kernel(self):
        length = self.L.m_as('m')
        k = self.K1.m_as('m**-2')
        if k == 0.0:
            self._kernel = drift5
            self._kernel_arguments = [_np.array([length])]
        else:
            if k > 0.0:
                k = _np.sqrt(k)
                kl = k * length
                s = _np.sin(kl)
                c = _np.cos(kl)
                sh = _np.sinh(kl)
                ch = _np.cosh(kl)
                self._kernel_arguments = [_np.array([
                    [c, (1 / k) * s, 0, 0, 0],
                    [-k * s, c, 0, 0, 0],
                    [0, 0, ch, (1 / k) * sh, 0],
                    [0, 0, k * sh, ch, 0],
                    [0, 0, 0, 0, 1]
                ])]
            else:
                k = _np.sqrt(-k)
                kl = k * length
                s = _np.sin(kl)
                c = _np.cos(kl)
                sh = _np.sinh(kl)
                ch = _np.cosh(kl)
                self._kernel_arguments = [_np.array([
                    [ch, (1 / k) * sh, 0, 0, 0],
                    [k * sh, ch, 0, 0, 0],
                    [0, 0, c, (1 / k) * s, 0],
                    [0, 0, -k * s, c, 0],
                    [0, 0, 0, 0, 1]
                ])]
            if self.integrator == FirstOrderIntegrator:
                self._kernel = batched_vector_matrix
            elif self.integrator == SecondOrderIntegrator:
                self._kernel = batched_vector_matrix_tensor
                if k > 0.0:
                    self._kernel_arguments.append(_np.zeros((3, 3, 3)))
                else:
                    self._kernel_arguments.append(_np.zeros((3, 3, 3)))


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

    def select_kernel(self, **kwargs) -> Callable:
        if self.ANGLE.magnitude == 0.0 and self.K1.magnitude == 0.0:
            return drift5
        else:
            return batched_vector_matrix

    def build_kernel_arguments(self):
        length = self.L.m_as('m')
        if self.ANGLE.magnitude == 0.0 and self.K1.magnitude == 0.0:
            return _np.array([length])
        elif self.ANGLE.magnitude == 0.0 and self.K1.magnitude != 0.0:
            return batched_vector_matrix
        else:
            theta = self.ANGLE.m_as('radian')
            k_bend = (theta/length)**2
            return batched_vector_matrix


class SBend(_ManzoniElement):
    pass


class RBend(_ManzoniElement):
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
