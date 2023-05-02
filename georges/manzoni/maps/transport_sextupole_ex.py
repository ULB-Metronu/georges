import numpy as np
from numba import njit
from numba.typed import List as nList


@njit(cache=True)
def compute_transport_sextupole_ex_matrix(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    R = np.zeros((6, 6))
    R[4, 4] = 1
    R[5, 5] = 1
    R[0, 0] = 1
    R[0, 1] = L
    R[1, 1] = 1
    R[2, 2] = 1
    R[2, 3] = L
    R[3, 3] = 1
    return R


@njit(cache=True)
def compute_transport_sextupole_ex_tensor(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    k2: float = element_parameters[1]
    d: float = element_parameters[len(element_parameters) - 1]
    k2 = k2 / (1 + d)
    T = np.zeros((6, 6, 6))
    T[0, 0, 0] = -(L**2) * k2 / 2
    T[0, 0, 1] = -(L**3) * k2 / 3
    T[0, 1, 1] = -(L**4) * k2 / 12
    T[0, 2, 2] = L**2 * k2 / 2
    T[0, 2, 3] = L**3 * k2 / 3
    T[0, 3, 3] = L**4 * k2 / 12
    T[1, 0, 0] = -L * k2
    T[1, 0, 1] = -(L**2) * k2
    T[1, 1, 1] = -(L**3) * k2 / 3
    T[1, 2, 2] = L * k2
    T[1, 2, 3] = L**2 * k2
    T[1, 3, 3] = L**3 * k2 / 3
    T[2, 0, 2] = L**2 * k2
    T[2, 0, 3] = L**3 * k2 / 3
    T[2, 1, 2] = L**3 * k2 / 3
    T[2, 1, 3] = L**4 * k2 / 6
    T[3, 0, 2] = 2 * L * k2
    T[3, 0, 3] = L**2 * k2
    T[3, 1, 2] = L**2 * k2
    T[3, 1, 3] = 2 * L**3 * k2 / 3
    return T
