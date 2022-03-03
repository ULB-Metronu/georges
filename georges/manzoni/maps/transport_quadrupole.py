from numba import njit
import numpy as np
from numba.typed import List as nList


@njit
def compute_transport_quadrupole_matrix(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    K1: float = element_parameters[1]
    R = np.zeros((6, 6))
    R[4, 4] = 1
    R[5, 5] = 1
    if K1 == 0:
        R[0, 0] = 1
        R[0, 1] = L
        R[1, 1] = 1
        R[2, 2] = 1
        R[2, 3] = L
        R[3, 3] = 1
    if K1 > 0:
        R[0, 0] = np.cos(np.sqrt(K1) * L)
        R[0, 1] = np.sin(np.sqrt(K1) * L) / np.sqrt(K1)
        R[1, 0] = -np.sqrt(K1) * np.sin(np.sqrt(K1) * L)
        R[1, 1] = np.cos(np.sqrt(K1) * L)
        R[2, 2] = np.cosh(np.sqrt(K1) * L)
        R[2, 3] = np.sinh(np.sqrt(K1) * L) / np.sqrt(K1)
        R[3, 2] = np.sqrt(K1) * np.sinh(np.sqrt(K1) * L)
        R[3, 3] = np.cosh(np.sqrt(K1) * L)
    if K1 < 0:
        R[0, 0] = np.cosh(L * np.sqrt(-K1))
        R[0, 1] = np.sinh(L * np.sqrt(-K1)) / np.sqrt(-K1)
        R[1, 0] = np.sqrt(-K1) * np.sinh(L * np.sqrt(-K1))
        R[1, 1] = np.cosh(L * np.sqrt(-K1))
        R[2, 2] = np.cos(L * np.sqrt(-K1))
        R[2, 3] = np.sin(L * np.sqrt(-K1)) / np.sqrt(-K1)
        R[3, 2] = -np.sqrt(-K1) * np.sin(L * np.sqrt(-K1))
        R[3, 3] = np.cos(L * np.sqrt(-K1))
    return R


@njit
def compute_transport_quadrupole_tensor(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    K1: float = element_parameters[1]
    T = np.zeros((6, 6, 6))
    if K1 > 0:
        T[0, 0, 4] = np.sqrt(K1) * L * np.sin(np.sqrt(K1) * L) / 2
        T[0, 1, 4] = -L * np.cos(np.sqrt(K1) * L) / 2 + np.sin(np.sqrt(K1) * L) / (2 * np.sqrt(K1))
        T[1, 0, 4] = np.sqrt(K1) * np.sin(np.sqrt(K1) * L) / 2 + K1 * L * np.cos(np.sqrt(K1) * L) / 2
        T[1, 1, 4] = np.sqrt(K1) * L * np.sin(np.sqrt(K1) * L) / 2
        T[2, 2, 4] = -np.sqrt(K1) * L * np.sinh(np.sqrt(K1) * L) / 2
        T[2, 3, 4] = -L * np.cosh(np.sqrt(K1) * L) / 2 + np.sinh(np.sqrt(K1) * L) / (2 * np.sqrt(K1))
        T[3, 2, 4] = -np.sqrt(K1) * np.sinh(np.sqrt(K1) * L) / 2 - K1 * L * np.cosh(np.sqrt(K1) * L) / 2
        T[3, 3, 4] = -np.sqrt(K1) * L * np.sinh(np.sqrt(K1) * L) / 2
    if K1 < 0:
        T[0, 0, 4] = -L * np.sqrt(-K1) * np.sinh(L * np.sqrt(-K1)) / 2
        T[0, 1, 4] = -(-L * np.sqrt(-K1) * np.cosh(L * np.sqrt(-K1)) + np.sinh(L * np.sqrt(-K1))) / (2 * np.sqrt(-K1))
        T[1, 0, 4] = K1 * L * np.cosh(L * np.sqrt(-K1)) / 2 - np.sqrt(-K1) * np.sinh(L * np.sqrt(-K1)) / 2
        T[1, 1, 4] = L * np.sqrt(-K1) * np.sinh(L * np.sqrt(-K1)) / 2
        T[2, 2, 4] = L * np.sqrt(-K1) * np.sin(L * np.sqrt(-K1)) / 2
        T[2, 3, 4] = L * np.cos(L * np.sqrt(-K1)) / 2 - np.sin(L * np.sqrt(-K1)) / (2 * np.sqrt(-K1))
        T[3, 2, 4] = -K1 * L * np.cos(L * np.sqrt(-K1)) / 2 + np.sqrt(-K1) * np.sin(L * np.sqrt(-K1)) / 2
        T[3, 3, 4] = -L * np.sqrt(-K1) * np.sin(L * np.sqrt(-K1)) / 2
    return T
