"""
TODO
"""
import numpy as np
from georges_core import Kinematics as _Kinematics


class Beam:
    def __init__(self, kinematics: _Kinematics, distribution: np.ndarray):
        self._kinematics = kinematics
        self._distribution = distribution.astype(float, copy=False)

    @classmethod
    def compute_pt(cls, dpp: np.ndarray, beta: float, first_order: bool = False) -> np.ndarray:
        """

        Args:
            dpp:
            beta: the relativistic beta
            first_order: approximate p_t from dpp at first order (valid only for beta close to 1)

        Returns:

        """
        if first_order:
            return dpp * beta
        else:
            return (-(2 / beta) + np.sqrt((2 / beta) ** 2 + 4 * (dpp**2 + 2 * dpp))) / 2

    @classmethod
    def compute_dpp(cls, pt: np.ndarray, beta: float, first_order: bool = False) -> np.ndarray:
        """

        Args:
            pt:
            beta: the relativistic beta
            first_order: approximate dpp from p_t at first order (valid only for beta close to 1)

        Returns:

        """
        if first_order:
            return pt / beta
        else:
            return np.sqrt(pt**2 + 2 * pt / beta + 1) - 1

    @property
    def kinematics(self) -> _Kinematics:
        return self._kinematics

    @property
    def distribution(self) -> np.ndarray:
        return self._distribution

    @distribution.setter
    def distribution(self, dist: np.ndarray):
        self._distribution = dist


class MadXBeam(Beam):
    def __init__(self, kinematics: _Kinematics, distribution: np.ndarray, first_order: bool = False):
        super().__init__(kinematics=kinematics, distribution=distribution)
        pt = Beam.compute_pt(
            distribution[:, -1],
            kinematics.beta,
            first_order=first_order,
        )
        self._distribution = np.zeros((len(distribution), 6))
        self._distribution[:, :-1] = distribution
        self._distribution[:, -1] = pt


class TransportBeam(Beam):
    def __init__(self, kinematics: _Kinematics, distribution: np.ndarray):
        super().__init__(kinematics=kinematics, distribution=distribution)
