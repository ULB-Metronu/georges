"""
TODO
"""
import numpy as np
from georges_core import Kinematics as _Kinematics


class Beam:
    def __init__(self, kinematics: _Kinematics, distribution: np.ndarray):
        self._kinematics = kinematics
        self._distribution = distribution

    @property
    def kinematics(self) -> _Kinematics:
        return self._kinematics

    @property
    def distribution(self) -> np.ndarray:
        return self._distribution
