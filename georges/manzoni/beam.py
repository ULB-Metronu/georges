"""
TODO
"""
import numpy as np
from georges_core import Kinematics as _Kinematics


class Beam:
    def __init__(self, kinematics: _Kinematics, distribution: np.ndarray):
        self._kinematics = kinematics
        self._distribution = distribution

    @classmethod
    def compute_pt(cls, dpp: np.ndarray, beta: float, first_order: bool = False) -> np.ndarray:
        """

        Args:
            dpp:
            beta:
            first_order:

        Returns:

        """
        if first_order:
            return dpp * beta
        else:
            return (- (2/beta) + np.sqrt((2/beta)**2 + 4 * (dpp**2 + 2 * dpp))) / 2

    @classmethod
    def compute_dpp(cls, pt: np.ndarray, beta: float, first_order: bool = False) -> np.ndarray:
        """

        Args:
            pt:
            beta:
            first_order:

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
