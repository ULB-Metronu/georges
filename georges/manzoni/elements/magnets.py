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
    """Define a Marker"""

    INTEGRATOR = None

    def propagate(self, beam_in: _np.ndarray, beam_out: _np.ndarray = None, parameters: list = None):
        _np.copyto(dst=beam_out, src=beam_in, casting="no")
        return beam_in, beam_out


class Matrix(_ManzoniElement):
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
    PARAMETERS = {
        "L": (0.0 * _ureg.m, "Drift length."),
        "APERTYPE": (None, "Aperture type (CIRCULAR, ELIPTIC or RECTANGULAR)"),
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
    PARAMETERS = {
        "APERTYPE": (None, "Aperture type (CIRCULAR, ELIPTIC or RECTANGULAR)"),
        "APERTURE": ([], ""),
    }
    """Parameters of the element, with their default value and their descriptions."""


class SRotation(_ManzoniElement):
    PARAMETERS = {
        "ANGLE": (0.0 * _ureg.radian, "Angle of rotation along the s-axis."),
    }
    """Parameters of the element, with their default value and their descriptions."""

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
        "APERTYPE": (None, "Aperture type (CIRCULAR, ELIPTIC or RECTANGULAR)"),
        "APERTURE": ([], ""),
        "KINEMATICS": (None, "Reference kinematics"),
    }


class Quadrupole(Magnet):
    PARAMETERS = {
        "L": (0.0 * _ureg.m, "Quadrupole length."),
        "K1": (0.0 * _ureg.m**-2, "Normalized gradient."),
        "K1S": (0.0 * _ureg.m**-2, "Normalized skew gradient."),
        "TILT": (0.0 * _ureg.radian, "Magnet tilt angle."),
    }
    """Parameters of the element, with their default value and their descriptions."""

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
    """Parameters of the element, with their default value and their descriptions."""

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
    pass


class RBend(Bend):
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
    PARAMETERS = {
        "H": (0.0 * _ureg.m**-1, "Inverse of the curvature radius."),
        "E1": (0.0 * _ureg.radian, "Entrance face angle."),
        "HGAP": (0.0 * _ureg.m, "Magnet gap."),
        "FINT": (0.0, "Fringe field integral."),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> list:
        h = self.H.m_as("m**-1")
        e1 = self.E1.m_as("radian")
        hgap = self.HGAP.m_as("m")
        fint = self.FINT
        return list(Bend.compute_fringe(h, e1, hgap, fint))


class Solenoid(Magnet):
    pass


class Multipole(_ManzoniElement):
    PARAMETERS = {
        "L": (0.0 * _ureg.m, "Magnet length."),
        "K1": (0.0 * _ureg.m**-2, "Quadrupolar normalized gradient."),
        "K2": (0.0 * _ureg.m**-3, "Sextupolar normalized gradient."),
    }
    """Parameters of the element, with their default value and their descriptions."""

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


class Sextupole(_ManzoniElement):
    PARAMETERS = {
        "L": (0.0 * _ureg.m, "Magnet length."),
        "K2": (0.0 * _ureg.m**-3, "Sextupolar normalized gradient."),
    }
    """Parameters of the element, with their default value and their descriptions."""

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


class Octupole(_ManzoniElement):
    pass


class Decapole(_ManzoniElement):
    pass


class Dodecapole(_ManzoniElement):
    pass


class Kicker(Magnet):
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
