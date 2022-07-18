from numba import njit
import numpy as np
from numba.typed import List as nList


@njit(parallel=True, fastmath=True)
def compute_transport_fringe_in_Ex_matrix(element_parameters: nList) -> np.ndarray:
    h: float = element_parameters[0]
    beta1: float = element_parameters[3]
    gap: float = element_parameters[5]
    fint: float = element_parameters[6]
    d: float = element_parameters[8]
    h = h / (1 + d)
    phi1 = fint * h * gap * np.cos(beta1) ** (-1) * (1 + np.sin(beta1) ** 2)
    R = np.zeros((6, 6))
    R[0, 0] = 1
    R[1, 0] = h * np.tan(beta1)
    R[1, 1] = 1
    R[2, 2] = 1
    R[3, 2] = -h * np.tan(beta1 - phi1)
    R[3, 3] = 1
    R[4, 4] = 1
    R[5, 5] = 1
    return R


@njit(parallel=True, fastmath=True)
def compute_transport_fringe_in_Ex_tensor(element_parameters: nList) -> np.ndarray:
    h: float = element_parameters[0]
    k1: float = element_parameters[1]
    beta1: float = element_parameters[3]
    gap: float = element_parameters[5]
    fint: float = element_parameters[6]
    r1: float = element_parameters[7]
    d: float = element_parameters[8]
    h = h / (1 + d)
    k1 = k1 / (1 + d)
    T = np.zeros((6, 6, 6))
    phi1 = fint * h * gap * np.cos(beta1) ** (-1) * (1 + np.sin(beta1) ** 2)
    T[0, 0, 0] = -h * np.tan(beta1) ** 2 / 2
    T[0, 2, 2] = h * np.cos(beta1) ** (-1) ** 2 / 2
    T[1, 0, 0] = k1 * np.tan(beta1) + (h * np.cos(beta1) ** (-1)) ** 3 / (2 * r1)
    T[1, 0, 1] = h * np.tan(beta1) ** 2
    T[1, 0, 4] = -h * np.tan(beta1)
    T[1, 2, 2] = -k1 * np.tan(beta1) + 0.5 * h ** 2 * np.tan(beta1) + h ** 2 + np.tan(beta1) ** 3 - h * np.cos(
        beta1) ** (-1) ** 3 / (2 * r1)
    T[1, 2, 3] = -h * np.tan(beta1) ** 2
    T[2, 0, 2] = h * np.tan(beta1) ** 2
    T[3, 0, 2] = -2 * k1 * np.tan(beta1) - h * np.cos(beta1) ** (-1) ** 3 / r1
    T[3, 0, 3] = -h * np.tan(beta1) ** 2
    T[3, 1, 2] = -h * np.cos(beta1) ** (-1) ** 2
    T[3, 2, 4] = -h * phi1 * np.cos(beta1 - phi1) ** (-2) + h * np.tan(beta1)
    return T
