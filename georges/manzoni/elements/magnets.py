"""
TODO
"""
from typing import Optional, Tuple

import numpy as _np

from ... import ureg as _ureg
from ..integrators import MadXIntegrator, MadXParaxialDriftIntegrator
from ..kernels import batched_vector_matrix
from .elements import ManzoniElement as _ManzoniElement


class Marker(_ManzoniElement):
    """
    Define a Marker
    """

    INTEGRATOR = None

    def propagate(self, beam_in: _np.ndarray, beam_out: _np.ndarray = None, parameters: list = None):
        _np.copyto(dst=beam_out, src=beam_in, casting="no")
        return beam_in, beam_out


class Matrix(_ManzoniElement):
    """
    Define a Matrix element.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Matrix with their default values.

    Examples:
        >>> m1 = Matrix('M1', MATRIX=_np.eye(6))
        >>> m1 #doctest: +NORMALIZE_WHITESPACE
            Matrix: {'NAME': 'M1',
            'AT_ENTRY': <Quantity(0, 'meter')>,
            'AT_CENTER': <Quantity(0, 'meter')>,
            'AT_EXIT': <Quantity(0, 'meter')>,
            'MATRIX': array([[1., 0., 0., 0., 0., 0.],
                             [0., 1., 0., 0., 0., 0.],
                             [0., 0., 1., 0., 0., 0.],
                             [0., 0., 0., 1., 0., 0.],
                             [0., 0., 0., 0., 1., 0.],
                             [0., 0., 0., 0., 0., 1.]])}
    """

    INTEGRATOR = None
    PARAMETERS = {
        "MATRIX": (_np.eye(6), "Transfer matrix."),
    }

    def propagate(
        self,
        beam_in: _np.ndarray,
        beam_out: Optional[_np.ndarray] = None,
        global_parameters: list = None,
    ) -> Tuple[_np.ndarray, _np.ndarray]:
        if self.integrator not in [MadXIntegrator, MadXParaxialDriftIntegrator]:
            beam_in[:, [-2, -1]] = beam_in[:, [-1, -2]]

        bvm = batched_vector_matrix(beam_in, beam_out, self.MATRIX)

        if self.integrator not in [MadXIntegrator, MadXParaxialDriftIntegrator]:
            bvm[0][:, [-1, -2]] = bvm[0][:, [-2, -1]]
            bvm[1][:, [-1, -2]] = bvm[1][:, [-2, -1]]

        return bvm

    @property
    def parameters(self) -> list:
        return [
            self.MATRIX,
        ]


class Gap(_ManzoniElement):
    """
    Define a gap where no physical geometry is placed

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Gap with their default values.

    Examples:
        >>> g1 = Gap('G1', L=1*_ureg.m)
        >>> g1 #doctest: +NORMALIZE_WHITESPACE
            Gap: {'NAME': 'G1',
            'AT_ENTRY': <Quantity(0, 'meter')>,
            'AT_CENTER': <Quantity(0, 'meter')>,
            'AT_EXIT': <Quantity(0, 'meter')>,
            'L': <Quantity(1, 'meter')>,
            'APERTYPE': None,
            'APERTURE': []}
    """

    PARAMETERS = {
        "L": (0.0 * _ureg.m, "Drift length."),
        "APERTYPE": (None, "Aperture type (CIRCULAR, ELLIPTICAL, RECTANGULAR or PHASE_SPACE)"),
        "APERTURE": ([], ""),
    }
    """Parameters of the element, with their default value and their descriptions."""

    def propagate(
        self,
        beam_in: _np.ndarray,
        beam_out: _np.ndarray = None,
        global_parameters: list = None,
    ) -> Tuple[_np.ndarray, _np.ndarray]:
        if self.L.magnitude == 0:
            _np.copyto(dst=beam_out, src=beam_in, casting="no")
            return beam_in, beam_out
        else:
            return self.integrator.propagate(self, beam_in, beam_out, global_parameters)

    @property
    def parameters(self) -> list:
        return list(
            map(
                float,
                [
                    self.L.m_as("meter"),
                ],
            ),
        )


class Drift(Gap):
    """
    Definition of a Drift.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Drift with their default values.

    Examples:
        >>> d1 = Drift('D1', L=1*_ureg.m, APERTYPE='CIRCULAR', APERTURE=[5*_ureg.cm])
        >>> d1 #doctest: +NORMALIZE_WHITESPACE
            Drift: {'NAME': 'D1',
                    'AT_ENTRY': <Quantity(0, 'meter')>,
                    'AT_CENTER': <Quantity(0, 'meter')>,
                    'AT_EXIT': <Quantity(0, 'meter')>,
                    'L': <Quantity(1, 'meter')>,
                    'APERTYPE': 'CIRCULAR',
                    'APERTURE': [<Quantity(5, 'centimeter')>]}
    """

    PARAMETERS = {
        "APERTYPE": (None, "Aperture type (CIRCULAR, ELLIPTICAL, RECTANGULAR or PHASE_SPACE)"),
        "APERTURE": ([], ""),
    }


class SRotation(_ManzoniElement):
    """
    Definition of a SRotation.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the SRotation with their default values.

    Examples:
        >>> s1 = SRotation('S1', ANGLE=10*_ureg.degrees)
        >>> s1 #doctest: +NORMALIZE_WHITESPACE
            SRotation: {'NAME': 'S1',
                        'AT_ENTRY': <Quantity(0, 'meter')>,
                        'AT_CENTER': <Quantity(0, 'meter')>,
                        'AT_EXIT': <Quantity(0, 'meter')>,
                        'ANGLE': <Quantity(10, 'degree')>}
    """

    PARAMETERS = {
        "ANGLE": (0.0 * _ureg.radian, "Angle of rotation along the s-axis."),
    }

    @property
    def parameters(self) -> list:
        return list(
            map(
                float,
                [
                    self.ANGLE.m_as("radian"),
                ],
            ),
        )


class Magnet(_ManzoniElement):
    PARAMETERS = {
        "APERTYPE": (None, "Aperture type (CIRCULAR, ELLIPTICAL, RECTANGULAR or PHASE_SPACE)"),
        "APERTURE": ([], ""),
        "KINEMATICS": (None, "Reference kinematics"),
    }


class Quadrupole(Magnet):
    """
    Define a Quadrupole magnet.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Quadrupole with their default values.

    Examples:
        >>> q1 = Quadrupole('Q1', L=1*_ureg.m, K1=3*_ureg.m**-2)
        >>> q1 #doctest: +NORMALIZE_WHITESPACE
            Quadrupole: {'NAME': 'Q1',
            'AT_ENTRY': <Quantity(0, 'meter')>,
            'AT_CENTER': <Quantity(0, 'meter')>,
            'AT_EXIT': <Quantity(0, 'meter')>,
            'APERTYPE': None, 'APERTURE': [],
            'KINEMATICS': None,
            'L': <Quantity(1, 'meter')>,
            'K1': <Quantity(3, '1 / meter ** 2')>,
            'K1S': <Quantity(0.0, '1 / meter ** 2')>,
            'TILT': <Quantity(0.0, 'radian')>}
    """

    PARAMETERS = {
        "L": (0.0 * _ureg.m, "Quadrupole length."),
        "K1": (0.0 * _ureg.m**-2, "Normalized gradient."),
        "K1S": (0.0 * _ureg.m**-2, "Normalized skew gradient."),
        "TILT": (0.0 * _ureg.radian, "Magnet tilt angle."),
    }

    @property
    def parameters(self) -> list:
        k1 = self.K1.m_as("m**-2")
        k1s = self.K1S.m_as("m**-2")
        tilt = self.TILT.m_as("radian")

        if k1s != 0.0 or tilt != 0.0:
            tilt -= _np.arctan2(k1s, k1) / 2
            k1 = _np.sqrt(k1**2 + k1s**2)

        return list(
            map(
                float,
                [
                    self.L.m_as("m"),
                    k1,
                    tilt,
                ],
            ),
        )


class Bend(Magnet):
    """
    Define a Bend magnet.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Bend with their default values.

    Examples:
        >>> b1 = Bend('B1', L=1*_ureg.m, ANGLE=30*_ureg.degrees, E1=5*_ureg.degrees, HGAP=2*_ureg.cm, K1=3*_ureg.m**-2)
        >>> b1 #doctest: +NORMALIZE_WHITESPACE
            Bend: {'NAME': 'B1',
            'AT_ENTRY': <Quantity(0, 'meter')>,
            'AT_CENTER': <Quantity(0, 'meter')>,
            'AT_EXIT': <Quantity(0, 'meter')>,
            'APERTYPE': None, 'APERTURE': [],
            'KINEMATICS': None,
            'ANGLE': <Quantity(30, 'degree')>,
            'K0': <Quantity(0.0, '1 / meter')>,
            'K1': <Quantity(3, '1 / meter ** 2')>,
            'K2': <Quantity(0.0, '1 / meter ** 3')>,
            'L': <Quantity(1, 'meter')>,
            'E1': <Quantity(5, 'degree')>,
            'E2': <Quantity(0.0, 'radian')>,
            'TILT': <Quantity(0.0, 'radian')>,
            'HGAP': <Quantity(2, 'centimeter')>,
            'FINT': 0.0,
            'FINTX': 0.0}
    """

    PARAMETERS = {
        "ANGLE": (0.0 * _ureg.radian, "Bending angle."),
        "K0": (0.0 * _ureg.m**-1, "Dipolar normalized gradient"),
        "K1": (0.0 * _ureg.m**-2, "Quadrupolar normalized gradient."),
        "K2": (0.0 * _ureg.m**-3, "Sextupolar normalized gradient."),
        "L": (0.0 * _ureg.m, "Magnet length."),
        "E1": (0.0 * _ureg.radian, "Entrance face angle."),
        "E2": (0.0 * _ureg.radian, "Exit face angle."),
        "TILT": (0.0 * _ureg.radian, "Magnet tilt angle."),
        "HGAP": (0.0 * _ureg.m, "Magnet gap."),
        "FINT": (0.0, "Fringe field integral."),
        "FINTX": (0.0, "Exit fringe field integral."),
    }

    @staticmethod
    def compute_fringe(h: float, e: float, hgap: float, fint: float) -> Tuple[float, float]:
        fringe_x = h * _np.tan(e)
        corr = h * (2 * hgap) * fint
        psi = e - corr * (1.0 / _np.cos(e)) * (1 + _np.sin(e) ** 2)
        fringe_y = -h * _np.tan(psi)
        return fringe_x, fringe_y

    @property
    def length(self) -> float:
        return self.L.m_as("m")

    @property
    def edges(self) -> Tuple[float, float]:
        return self.E1.m_as("radian"), self.E2.m_as("radian")

    @property
    def fringe_field_integrals(self) -> Tuple[float, float]:
        fint = self.FINT
        fintx = self.FINTX if self.FINTX >= 0 else fint  # For exact compatibility with MAD-X
        return fint, fintx

    @property
    def parameters(self) -> list:
        # Generic parameters
        length = self.length
        h = self.ANGLE.m_as("radian") / self.length
        if self.K0.m_as("m**-1") is None or self.K0.m_as("m**-1") == 0:
            k0 = h
        else:
            k0 = self.K0.m_as("m**-1")
        hgap = self.HGAP.m_as("m")
        e1, e2 = self.edges
        fint, fintx = self.fringe_field_integrals
        entrance_fringe_x, entrance_fringe_y = Bend.compute_fringe(h, e1, hgap, fint)
        exit_fringe_x, exit_fringe_y = Bend.compute_fringe(h, e2, hgap, fintx)

        return list(
            map(
                float,
                [
                    length,  # 0
                    self.ANGLE.m_as("radian"),  # 1
                    self.K1.m_as("m**-2"),  # 2
                    self.K2.m_as("m**-3"),  # 3
                    self.TILT.m_as("radian"),  # 4
                    h,  # 5
                    k0,  # 6
                    entrance_fringe_x,  # 7
                    entrance_fringe_y,  # 8
                    exit_fringe_x,  # 9
                    exit_fringe_y,  # 10
                ],
            ),
        )


class SBend(Bend):
    """
    Definition of a SBend
    """

    pass


class RBend(Bend):
    """
    Definition of a RBend
    """

    @property
    def length(self) -> float:
        length = self.L.m_as("m")
        angle = self.ANGLE.m_as("rad")
        if angle > 1e-8:
            return length * angle / (2.0 * _np.sin(angle / 2.0))
        else:
            return length

    @property
    def edges(self) -> Tuple[float, float]:
        alpha = self.ANGLE.m_as("radian") / 2.0
        return self.E1.m_as("radian") + alpha, self.E2.m_as("radian") + alpha


class Fringein(Bend):
    """
    Define a Fringe field at the entrance of a bending magnet.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Fringein with their default values.

    Examples:
        >>> f1 = Fringein('F1', L=1*_ureg.mm, ANGLE=30*_ureg.degrees, E1=5*_ureg.degrees,
        ...               HGAP=2*_ureg.cm, K1=3*_ureg.m**-2, R1=1*_ureg.cm)
        >>> f1 #doctest: +NORMALIZE_WHITESPACE
            Fringein: {'NAME': 'F1',
                       'AT_ENTRY': <Quantity(0, 'meter')>,
                       'AT_CENTER': <Quantity(0, 'meter')>,
                       'AT_EXIT': <Quantity(0, 'meter')>,
                       'APERTYPE': None, 'APERTURE': [],
                       'KINEMATICS': None,
                       'ANGLE': <Quantity(30, 'degree')>,
                       'K0': <Quantity(0.0, '1 / meter')>,
                       'K1': <Quantity(3, '1 / meter ** 2')>,
                       'K2': <Quantity(0.0, '1 / meter ** 3')>,
                       'L': <Quantity(1, 'millimeter')>,
                       'E1': <Quantity(5, 'degree')>,
                       'E2': <Quantity(0.0, 'radian')>,
                       'TILT': <Quantity(0.0, 'radian')>,
                       'HGAP': <Quantity(2, 'centimeter')>,
                       'FINT': 0.0, 'FINTX': 0.0,
                       'R1': <Quantity(1, 'centimeter')>}
    """

    PARAMETERS = {
        "ANGLE": (0.0 * _ureg.radian, "Bending angle."),
        "L": (0.0 * _ureg.m, "Magnet length."),
        "K1": (0.0 * _ureg.m**-2, "Quadrupolar normalized gradient."),
        "E1": (0.0 * _ureg.radian, "Entrance face angle."),
        "HGAP": (0.0 * _ureg.m, "Magnet gap."),
        "FINT": (0.0, "Fringe field integral."),
        "R1": (_np.inf * _ureg.m, "Entrance face curvature radius"),
    }

    @property
    def parameters(self) -> list:
        # Generic parameters
        return list(
            map(
                float,
                [
                    self.length,  # 0
                    self.ANGLE.m_as("radian"),  # 1
                    self.K1.m_as("m**-2"),  # 2
                    self.E1.m_as("radian"),  # 3
                    self.HGAP.m_as("m"),  # 4
                    self.FINT,  # 5
                    self.R1.m_as("m"),  # 6
                ],
            ),
        )


class Fringeout(Bend):
    """
    Define a Fringe field at the exit of a bending magnet.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Fringeout with their default values.

    Examples:
        >>> f2 = Fringeout('F2', L=1*_ureg.mm, ANGLE=30*_ureg.degrees, E1=5*_ureg.degrees,
        ...                 HGAP=2*_ureg.cm, K1=3*_ureg.m**-2, R2=1*_ureg.cm)
        >>> f2 #doctest: +NORMALIZE_WHITESPACE
            Fringeout: {'NAME': 'F2',
                       'AT_ENTRY': <Quantity(0, 'meter')>,
                       'AT_CENTER': <Quantity(0, 'meter')>,
                       'AT_EXIT': <Quantity(0, 'meter')>,
                       'APERTYPE': None, 'APERTURE': [],
                       'KINEMATICS': None,
                       'ANGLE': <Quantity(30, 'degree')>,
                       'K0': <Quantity(0.0, '1 / meter')>,
                       'K1': <Quantity(3, '1 / meter ** 2')>,
                       'K2': <Quantity(0.0, '1 / meter ** 3')>,
                       'L': <Quantity(1, 'millimeter')>,
                       'E1': <Quantity(5, 'degree')>,
                       'E2': <Quantity(0.0, 'radian')>,
                       'TILT': <Quantity(0.0, 'radian')>,
                       'HGAP': <Quantity(2, 'centimeter')>,
                       'FINT': 0.0, 'FINTX': 0.0,
                       'R2': <Quantity(1, 'centimeter')>}
    """

    PARAMETERS = {
        "ANGLE": (0.0 * _ureg.radian, "Bending angle."),
        "L": (0.0 * _ureg.m, "Magnet length."),
        "K1": (0.0 * _ureg.m**-2, "Quadrupolar normalized gradient."),
        "E2": (0.0 * _ureg.radian, "Exit face angle."),
        "HGAP": (0.0 * _ureg.m, "Magnet gap."),
        "FINTX": (0.0, "Exit fringe field integral."),
        "R2": (_np.inf * _ureg.m, "Entrance face curvature radius"),
    }

    @property
    def parameters(self) -> list:
        return list(
            map(
                float,
                [
                    self.length,  # 0
                    self.ANGLE.m_as("radian"),  # 1
                    self.K1.m_as("m**-2"),  # 2
                    self.E2.m_as("radian"),  # 3
                    self.HGAP.m_as("m"),  # 4
                    self.FINTX,  # 5
                    self.R2.m_as("m"),  # 6
                ],
            ),
        )


class DipEdge(Magnet):
    """
    Define a DipEdge element.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the DipEdge with their default values.

    Examples:
        >>> d1 = DipEdge('D1', H=1*_ureg.mm**-1)
        >>> d1 #doctest: +NORMALIZE_WHITESPACE
            DipEdge: {'NAME': 'D1',
                      'AT_ENTRY': <Quantity(0, 'meter')>,
                      'AT_CENTER': <Quantity(0, 'meter')>,
                      'AT_EXIT': <Quantity(0, 'meter')>,
                      'APERTYPE': None,
                      'APERTURE': [],
                      'KINEMATICS': None,
                      'H': <Quantity(1, '1 / millimeter')>,
                      'E1': <Quantity(0.0, 'radian')>,
                      'HGAP': <Quantity(0.0, 'meter')>,
                      'FINT': 0.0}
    """

    PARAMETERS = {
        "H": (0.0 * _ureg.m**-1, "Inverse of the curvature radius."),
        "E1": (0.0 * _ureg.radian, "Entrance face angle."),
        "HGAP": (0.0 * _ureg.m, "Magnet gap."),
        "FINT": (0.0, "Fringe field integral."),
    }

    @property
    def parameters(self) -> list:
        h = self.H.m_as("m**-1")
        e1 = self.E1.m_as("radian")
        hgap = self.HGAP.m_as("m")
        fint = self.FINT
        return list(Bend.compute_fringe(h, e1, hgap, fint))


class Solenoid(Magnet):
    pass


class Multipole(Magnet):
    """
    Define a Multipole magnet.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Multipole with their default values.

    Examples:
        >>> m1 = Multipole('M1', L=10*_ureg.cm, K1=4*_ureg.m**-2, K2=0.1*_ureg.m**-3)
        >>> m1 #doctest: +NORMALIZE_WHITESPACE
            Multipole: {'NAME': 'M1',
                        'AT_ENTRY': <Quantity(0, 'meter')>,
                        'AT_CENTER': <Quantity(0, 'meter')>,
                        'AT_EXIT': <Quantity(0, 'meter')>,
                        'APERTYPE': None,
                        'APERTURE': [],
                        'KINEMATICS': None,
                        'L': <Quantity(10, 'centimeter')>,
                        'K1': <Quantity(4, '1 / meter ** 2')>,
                        'K2': <Quantity(0.1, '1 / meter ** 3')>}
    """

    PARAMETERS = {
        "L": (0.0 * _ureg.m, "Magnet length."),
        "K1": (0.0 * _ureg.m**-2, "Quadrupolar normalized gradient."),
        "K2": (0.0 * _ureg.m**-3, "Sextupolar normalized gradient."),
    }

    @property
    def parameters(self) -> list:
        return list(
            map(
                float,
                [
                    self.L.m_as("m"),
                    self.K1.m_as("m**-2"),
                    self.K2.m_as("m**-3"),
                ],
            ),
        )


class Sextupole(Magnet):
    """
    Define a Sextupole magnet.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Sextupole with their default values.

    Examples:
        >>> s1 = Sextupole('S1', L=1*_ureg.m, K2=3*_ureg.m**-3)
        >>> s1 #doctest: +NORMALIZE_WHITESPACE
            Sextupole: {'NAME': 'S1',
            'AT_ENTRY': <Quantity(0, 'meter')>,
            'AT_CENTER': <Quantity(0, 'meter')>,
            'AT_EXIT': <Quantity(0, 'meter')>,
            'APERTYPE': None, 'APERTURE': [],
            'KINEMATICS': None,
            'L': <Quantity(1, 'meter')>,
            'K2': <Quantity(3, '1 / meter ** 3')>}
    """

    PARAMETERS = {
        "L": (0.0 * _ureg.m, "Magnet length."),
        "K2": (0.0 * _ureg.m**-3, "Sextupolar normalized gradient."),
    }

    @property
    def parameters(self) -> list:
        return list(
            map(
                float,
                [
                    self.L.m_as("m"),
                    self.K2.m_as("m**-3"),
                ],
            ),
        )


class Octupole(Magnet):
    pass


class Decapole(Magnet):
    pass


class Dodecapole(Magnet):
    pass


class Kicker(Magnet):
    """
    Define a Kicker element.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Kicker with their default values.

    Examples:
        >>> k1 = Kicker('K1', L=10*_ureg.cm, HKICK=0.1, VKICK=-0.2)
        >>> k1 #doctest: +NORMALIZE_WHITESPACE
            Kicker: {'NAME': 'K1',
                     'AT_ENTRY': <Quantity(0, 'meter')>,
                     'AT_CENTER': <Quantity(0, 'meter')>,
                     'AT_EXIT': <Quantity(0, 'meter')>,
                     'APERTYPE': None, 'APERTURE': [],
                     'KINEMATICS': None,
                     'L': <Quantity(10, 'centimeter')>,
                     'HKICK': 0.1,
                     'VKICK': -0.2,
                     'TILT': <Quantity(0.0, 'radian')>}
    """

    PARAMETERS = {
        "L": (0.0 * _ureg.m, "Kicker length."),
        "HKICK": (0.0, "The momentum change in the horizontal plane."),
        "VKICK": (0.0, "The momentum change in the vertical plane."),
        "TILT": (
            0.0 * _ureg.radian,
            "The roll angle about the longitudinal axis. A positive angle represents a "
            "clockwise rotation of the kicker",
        ),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> list:
        tilt = -self.TILT.m_as("radian")
        ct = _np.cos(tilt)
        st = _np.sin(tilt)
        hkick = self.HKICK
        vkick = self.VKICK
        _hkick = ct * hkick + st * vkick
        _vkick = ct * vkick - st * hkick
        return list(
            map(
                float,
                [
                    self.L.m_as("m"),
                    _hkick,
                    _vkick,
                ],
            ),
        )


class TKicker(Kicker):
    pass


class HKicker(Magnet):
    """
    Define a HKicker element.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the HKicker with their default values.

    Examples:
        >>> h1 = HKicker('H1', L=10*_ureg.cm, KICK=0.1)
        >>> h1 #doctest: +NORMALIZE_WHITESPACE
            HKicker: {'NAME': 'H1',
                     'AT_ENTRY': <Quantity(0, 'meter')>,
                     'AT_CENTER': <Quantity(0, 'meter')>,
                     'AT_EXIT': <Quantity(0, 'meter')>,
                     'APERTYPE': None, 'APERTURE': [],
                     'KINEMATICS': None,
                     'L': <Quantity(10, 'centimeter')>,
                     'KICK': 0.1,
                     'TILT': <Quantity(0.0, 'radian')>}
    """

    PARAMETERS = {
        "L": (0.0 * _ureg.m, "Kicker length."),
        "KICK": (0.0, "The momentum change."),
        "TILT": (
            0.0 * _ureg.radian,
            "The roll angle about the longitudinal axis. A positive angle represents a "
            "clockwise rotation of the kicker",
        ),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> list:
        tilt = -self.TILT.m_as("radian")
        ct = _np.cos(tilt)
        st = _np.sin(tilt)
        kick = self.KICK
        hkick = ct * kick
        vkick = -st * kick
        return list(
            map(
                float,
                [
                    self.L.m_as("m"),
                    hkick,
                    vkick,
                ],
            ),
        )


class VKicker(Magnet):
    """
    Define a VKicker element.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the VKicker with their default values.

    Examples:
        >>> v1 = VKicker('V1', L=10*_ureg.cm, KICK=-0.2)
        >>> v1 #doctest: +NORMALIZE_WHITESPACE
            VKicker: {'NAME': 'V1',
                     'AT_ENTRY': <Quantity(0, 'meter')>,
                     'AT_CENTER': <Quantity(0, 'meter')>,
                     'AT_EXIT': <Quantity(0, 'meter')>,
                     'APERTYPE': None, 'APERTURE': [],
                     'KINEMATICS': None,
                     'L': <Quantity(10, 'centimeter')>,
                     'KICK': -0.2,
                     'TILT': <Quantity(0.0, 'radian')>}
    """

    PARAMETERS = {
        "L": (0.0 * _ureg.m, "Kicker length."),
        "KICK": (0.0 * _ureg.radian, "The momentum change."),
        "TILT": (
            0.0 * _ureg.radian,
            "The roll angle about the longitudinal axis. A positive angle represents a "
            "clockwise rotation of the kicker",
        ),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> list:
        tilt = -self.TILT.m_as("radian")
        ct = _np.cos(tilt)
        st = _np.sin(tilt)
        kick = self.KICK
        hkick = st * kick
        vkick = ct * kick
        return list(
            map(
                float,
                [
                    self.L.m_as("m"),
                    hkick,
                    vkick,
                ],
            ),
        )
