import numpy as np
from numba import njit
from numba.typed import List as nList
from numpy import cos, cosh, sin, sinh, sqrt


@njit(cache=True)
def compute_transport_quadrupole_ex_matrix(element_parameters: nList) -> np.ndarray:

    L: float = element_parameters[0]
    k1: float = element_parameters[1]
    d: float = element_parameters[2]
    k1 = k1 / (1 + d)
    R = np.zeros((6, 6))
    R[4, 4] = 1
    R[5, 5] = 1
    if k1 == 0:
        R[0, 0] = 1
        R[0, 1] = L
        R[1, 1] = 1
        R[2, 2] = 1
        R[2, 3] = L
        R[3, 3] = 1
    if k1 > 0:
        R[0, 0] = cos(L * sqrt(k1))
        R[0, 1] = sin(L * sqrt(k1)) / sqrt(k1)
        R[1, 0] = -sqrt(k1) * sin(L * sqrt(k1))
        R[1, 1] = cos(L * sqrt(k1))
        R[2, 2] = cosh(L * sqrt(k1))
        R[2, 3] = sinh(L * sqrt(k1)) / sqrt(k1)
        R[3, 2] = sqrt(k1) * sinh(L * sqrt(k1))
        R[3, 3] = cosh(L * sqrt(k1))
    if k1 < 0:
        R[0, 0] = cosh(L * sqrt(-k1))
        R[0, 1] = sinh(L * sqrt(-k1)) / sqrt(-k1)
        R[1, 0] = sqrt(-k1) * sinh(L * sqrt(-k1))
        R[1, 1] = cosh(L * sqrt(-k1))
        R[2, 2] = cos(L * sqrt(-k1))
        R[2, 3] = sin(L * sqrt(-k1)) / sqrt(-k1)
        R[3, 2] = -sqrt(-k1) * sin(L * sqrt(-k1))
        R[3, 3] = cos(L * sqrt(-k1))
    return R


@njit(cache=True)
def compute_transport_quadrupole_ex_tensor(element_parameters: nList) -> np.ndarray:

    L: float = element_parameters[0]
    k1: float = element_parameters[1]
    d: float = element_parameters[len(element_parameters) - 1]
    k1 = k1 / (1 + d)
    T = np.zeros((6, 6, 6))
    if k1 > 0:
        T[0, 0, 5] = L * sqrt(k1) * sin(L * sqrt(k1)) / 2
        T[0, 1, 5] = -L * cos(L * sqrt(k1)) / 2 + sin(L * sqrt(k1)) / (2 * sqrt(k1))
        T[1, 0, 5] = L * k1 * cos(L * sqrt(k1)) / 2 + sqrt(k1) * sin(L * sqrt(k1)) / 2
        T[1, 1, 5] = L * sqrt(k1) * sin(L * sqrt(k1)) / 2
        T[2, 2, 5] = -L * sqrt(k1) * sinh(L * sqrt(k1)) / 2
        T[2, 3, 5] = -L * cosh(L * sqrt(k1)) / 2 + sinh(L * sqrt(k1)) / (2 * sqrt(k1))
        T[3, 2, 5] = -L * k1 * cosh(L * sqrt(k1)) / 2 - sqrt(k1) * sinh(L * sqrt(k1)) / 2
        T[3, 3, 5] = -L * sqrt(k1) * sinh(L * sqrt(k1)) / 2
    if k1 < 0:
        T[0, 0, 5] = -L * sqrt(-k1) * sinh(L * sqrt(-k1)) / 2
        T[0, 1, 5] = -L * cosh(L * sqrt(-k1)) / 2 + sinh(L * sqrt(-k1)) / (2 * sqrt(-k1))
        T[1, 0, 5] = L * k1 * cosh(L * sqrt(-k1)) / 2 - sqrt(-k1) * sinh(L * sqrt(-k1)) / 2
        T[1, 1, 5] = -L * sqrt(-k1) * sinh(L * sqrt(-k1)) / 2
        T[2, 2, 5] = L * sqrt(-k1) * sin(L * sqrt(-k1)) / 2
        T[2, 3, 5] = -L * cos(L * sqrt(-k1)) / 2 + sin(L * sqrt(-k1)) / (2 * sqrt(-k1))
        T[3, 2, 5] = -L * k1 * cos(L * sqrt(-k1)) / 2 + sqrt(-k1) * sin(L * sqrt(-k1)) / 2
        T[3, 3, 5] = L * sqrt(-k1) * sin(L * sqrt(-k1)) / 2
    return T
