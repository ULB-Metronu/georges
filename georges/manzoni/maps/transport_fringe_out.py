from numba import njit
import numpy as np
from numpy import cos, sin, tan
from numba.typed import List as nList


@njit(cache=True)
def compute_transport_fringe_out_matrix(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    beta2: float = element_parameters[3]
    gap: float = element_parameters[4]
    fintx: float = element_parameters[5]
    h = alpha / L
    R = np.zeros((6, 6))
    phi2 = fintx * h * gap * cos(beta2) ** (-1) * (1 + sin(beta2) ** 2)
    R[0, 0] = 1
    R[1, 0] = h * tan(beta2)
    R[1, 1] = 1
    R[2, 2] = 1
    R[3, 2] = -h * tan(beta2 - phi2)
    R[3, 3] = 1
    R[4, 4] = 1
    R[5, 5] = 1
    return R


@njit(cache=True)
def compute_transport_fringe_out_tensor(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    k1: float = element_parameters[2]
    beta2: float = element_parameters[3]
    gap: float = element_parameters[4]
    fintx: float = element_parameters[5]
    r2: float = element_parameters[6]
    h = alpha / L
    T = np.zeros((6, 6, 6))
    phi2 = fintx * h * gap * cos(beta2) ** (-1) * (1 + sin(beta2) ** 2)
    T[0, 0, 0] = h * tan(beta2) ** 2 / 2
    T[0, 2, 2] = -h * cos(beta2) ** (-1) ** 2 / 2
    T[1, 0, 0] = k1 * tan(beta2) - 0.5 * h ** 2 * tan(beta2) ** 3 + h * cos(beta2) ** (-1) ** 3 / (2 * r2)
    T[1, 0, 1] = -h * tan(beta2) ** 2
    T[1, 0, 4] = -h * tan(beta2)
    T[1, 2, 2] = -k1 * tan(beta2) - 0.5 * h ** 2 * tan(beta2) ** 3 - h * cos(beta2) ** (-1) ** 3 / (2 * r2)
    T[1, 2, 3] = h * tan(beta2) ** 2
    T[2, 0, 2] = -h * tan(beta2) ** 2
    T[3, 0, 2] = -2 * k1 * tan(beta2) + h ** 2 * tan(beta2) * cos(beta2) ** (-1) ** 2 - h * cos(beta2) ** (
        -1) ** 3 / r2
    T[3, 0, 3] = h * tan(beta2) ** 2
    T[3, 1, 2] = h * cos(beta2) ** (-1) ** 2
    T[3, 2, 4] = -h * phi2 * cos(beta2 - phi2) ** (-2) + h * tan(beta2)
    return T
