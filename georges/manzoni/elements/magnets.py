"""
TODO
"""
from typing import Optional, Callable
import numpy as _np
from ... import ureg as _ureg
from .elements import ManzoniElement as _ManzoniElement
from ..kernels import drift5, batched_vector_matrix


class Marker(_ManzoniElement):
    def propagate(self, beam: _np.ndarray, out: Optional[_np.ndarray] = None):
        _np.copyto(out, beam, 'no')
        return out


class Drift(_ManzoniElement):
    PARAMETERS = {
        'L': (0.0 * _ureg.m, 'Drift length.'),
    }
    """Parameters of the element, with their default value and their description ."""

    def build_kernel(self):
        self._kernel = drift5
        self._kernel_arguments = _np.array([self.L.m_as('m')])


class Rotation(_ManzoniElement):
    PARAMETERS = {
        'ANGLE': (0.0 * _ureg.radian, 'Angle of rotation along the s-axis.')
    }
    """Parameters of the element, with their default value and their description ."""

    def build_kernel(self):
        self._kernel = batched_vector_matrix
        angle = self.ANGLE.m_as('radian')
        self._kernel_arguments = _np.array([
            [_np.cos(angle), 0, -_np.sin(angle), 0, 0],
            [0, _np.cos(angle), 0, -_np.sin(angle), 0],
            [_np.sin(angle), 0, _np.cos(angle), 0, 0],
            [0, _np.sin(angle), 0, _np.cos(angle), 0],
            [0, 0, 0, 0, 1],
        ])


class Quadrupole(_ManzoniElement):
    PARAMETERS = {
        'L': (0.0 * _ureg.m, 'Quadrupole length.'),
        'K1': (0.0 * _ureg.m**-2, 'Normalized gradient.'),
        'K1L': (0.0 * _ureg.m**-1, 'Normalized integrated gradient.')
    }
    """Parameters of the element, with their default value and their description ."""

    def build_kernel(self):
        length = self.L.m_as('m')
        k = self.K1.m_as('m**-2')
        if k == 0.0:
            self._kernel = drift5
            self._kernel_arguments = _np.array([length])
        else:
            self._kernel = batched_vector_matrix
            if k > 0.0:
                k = _np.sqrt(k)
                kl = k * length
                s = _np.sin(kl)
                c = _np.cos(kl)
                sh = _np.sinh(kl)
                ch = _np.cosh(kl)
                self._kernel_arguments = _np.array([
                    [c, (1 / k) * s, 0, 0, 0],
                    [-k * s, c, 0, 0, 0],
                    [0, 0, ch, (1 / k) * sh, 0],
                    [0, 0, k * sh, ch, 0],
                    [0, 0, 0, 0, 1]
                ])
            else:
                k = _np.sqrt(-k)
                kl = k * length
                s = _np.sin(kl)
                c = _np.cos(kl)
                sh = _np.sinh(kl)
                ch = _np.cosh(kl)
                self._kernel_arguments = _np.array([
                    [ch, (1 / k) * sh, 0, 0, 0],
                    [k * sh, ch, 0, 0, 0],
                    [0, 0, c, (1 / k) * s, 0],
                    [0, 0, -k * s, c, 0],
                    [0, 0, 0, 0, 1]
                ])


class Bend(_ManzoniElement):
    PARAMETERS = {
        'ANGLE': (0.0 * _ureg.radian, 'Bending angle.'),
        'K1': (0.0 * _ureg.m**-2, 'Quadrupolar normalized gradient.'),
        'L': (0.0 * _ureg.m, 'Magnet length.'),
        'E1': (0.0 * _ureg.radian, 'Entrance face angle.'),
        'E2': (0.0 * _ureg.radian, 'Exirt face angle.'),
    }
    """Parameters of the element, with their default value and their description ."""

    def select_kernel(self, **kwargs) -> Callable:
        if self.ANGLE.magnitude == 0.0 and self.K1.magnitude == 0.0:
            return drift4
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




    # Definition of the main variables
    theta = e[INDEX_ANGLE]
    length = e[INDEX_LENGTH]
    e1 = e[INDEX_E1] + e1
    e2 = e[INDEX_E2] + e2
    k_bend = (theta/length)**2
    k = np.sqrt(np.abs(k_bend + e[INDEX_K1] / e[INDEX_BRHO]))
    k1 = np.sqrt(np.abs(e[INDEX_K1] / e[INDEX_BRHO]))
    kl = k * length
    k1l = k1 * length
    s_dpp = np.sin(theta)
    c_dpp = np.cos(theta)

    # Horizontal plane
    if k_bend + e[INDEX_K1] / e[INDEX_BRHO] > 0:
        s = np.sin(kl)
        c = np.cos(kl)
        m1 = [c, (1 / k) * s, 0, 0, (length / theta) * (1 - c_dpp)]
        m2 = [-k * s, c, 0, 0, s_dpp]
    elif k_bend + e[INDEX_K1] / e[INDEX_BRHO] < 0:
        sh = np.sinh(kl)
        ch = np.cosh(kl)
        m1 = [ch, (1 / k) * sh, 0, 0, (length / theta) * (1 - c_dpp)]
        m2 = [k * sh, ch, 0, 0, s_dpp]
    else:  # k_bend + e[INDEX_K1] == 0
        m1 = [1, length, 0, 0, (length / theta) * (1 - c_dpp)]
        m2 = [0, 1, 0, 0, s_dpp]

    # Vertical plane
    if e[INDEX_K1] < 0:
        s = np.sin(k1l)
        c = np.cos(k1l)
        m3 = [0, 0, c, (1 / k1) * s, 0]
        m4 = [0, 0, -k1 * s, c, 0]
    elif e[INDEX_K1] > 0:
        sh = np.sinh(k1l)
        ch = np.cosh(k1l)
        m3 = [0, 0, ch, (1 / k1) * sh, 0]
        m4 = [0, 0, k1 * sh, ch, 0]
    else:  # k1 == 0
        m3 = [0, 0, 1, length, 0]
        m4 = [0, 0, 0, 1, 0]

    # Construct 5D matrix
    m_b = np.stack([
        m1,
        m2,
        m3,
        m4,
        [0, 0, 0, 0, 1],
    ])

    # Poleface angle
    if e1 == 0 and e2 == 0:
        return m_b
    else:
        h = theta / length
        psi1 = e[INDEX_FINT] * h * e[INDEX_HGAP] * (1/np.cos(e1)) * (1 + (np.sin(e1))**2)
        psi2 = e[INDEX_FINT] * h * e[INDEX_HGAP] * (1/np.cos(e2)) * (1 + (np.sin(e2))**2)
        k1_x = h * np.tan(e1)
        k1_y = h * np.tan(e1-psi1)
        k2_x = h * np.tan(e2)
        k2_y = h * np.tan(e2-psi2)

        m_e1 = np.array(
            [
                [1, 0, 0, 0, 0],
                [k1_x, 1, 0, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, -k1_y, 1, 0],
                [0, 0, 0, 0, 1]
            ]
        )
        m_e2 = np.array(
            [
                [1, 0, 0, 0, 0],
                [k2_x, 1, 0, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, -k2_y, 1, 0],
                [0, 0, 0, 0, 1]
            ]
        )
        if multiply:
            return m_e2 @ m_b @ m_e1
        else:
            return [m_e2, m_b, m_e1]



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
