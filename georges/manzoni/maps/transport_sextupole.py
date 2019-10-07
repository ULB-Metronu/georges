from numba import njit
import numpy as np


@njit
def compute_transport_sextupole_matrix(L: float, *args) -> np.ndarray:

    R = np.zeros((6, 6))

    R[0, 0] = 1
    R[0, 1] = L
    R[1, 0] = 0
    R[1, 1] = 1
    R[2, 2] = 1
    R[2, 3] = L
    R[3, 2] = 0
    R[3, 3] = 1
    R[4, 4] = 1
    R[5, 5] = 1

    return R


@njit
def compute_transport_sextupole_tensor(L: float, K2: float, *args) -> np.ndarray:
    T = np.zeros((6, 6, 6))

    T[0, 0, 0] = -(1 / 2) * K2 * L ** 2
    T[0, 0, 1] = -(1 / 3) * K2 * L ** 3
    T[0, 1, 1] = -(1 / 12) * K2 * L ** 4
    T[0, 2, 2] = (1 / 2) * K2 * L ** 2
    T[0, 2, 3] = (1 / 3) * K2 * L ** 3
    T[0, 3, 3] = (1 / 12) * K2 * L ** 4
    T[1, 0, 0] = -K2 * L
    T[1, 0, 1] = -K2 * L ** 2
    T[1, 1, 1] = -(1 / 3) * K2 * L ** 3
    T[1, 2, 2] = K2 * L
    T[1, 2, 3] = K2 * L ** 2
    T[1, 3, 3] = (1 / 3) * K2 * L ** 3
    T[2, 0, 2] = K2 * L ** 2
    T[2, 0, 3] = (1 / 3) * K2 * L ** 3
    T[2, 1, 2] = (1 / 3) * K2 * L ** 3
    T[2, 1, 3] = (1 / 6) * K2 * L ** 4
    T[3, 0, 2] = 2 * K2 * L
    T[3, 0, 3] = K2 * L ** 2
    T[3, 1, 2] = K2 * L ** 2
    T[3, 1, 3] = (2 / 3) * K2 * L ** 3

    return T
