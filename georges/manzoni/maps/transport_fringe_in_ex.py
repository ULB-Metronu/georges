from numba import njit
import numpy as np
from numpy import cos, sin, tan
from numba.typed import List as nList


@njit(cache=True)
def compute_transport_fringe_in_ex_matrix(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    beta1: float = element_parameters[3]
    gap: float = element_parameters[4]
    fint: float = element_parameters[5]
    d: float = element_parameters[len(element_parameters) - 1]

    h = alpha / L
    h = h / (1 + d)
    phi1 = fint * h * gap * cos(beta1) ** (-1) * (1 + sin(beta1) ** 2)
    R = np.zeros((6, 6))
    R[0, 0] = 1
    R[1, 0] = h * tan(beta1)
    R[1, 1] = 1
    R[2, 2] = 1
    R[3, 2] = -h * tan(beta1 - phi1)
    R[3, 3] = 1
    R[4, 4] = 1
    R[5, 5] = 1
    return R


@njit(cache=True)
def compute_transport_fringe_in_ex_tensor(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    k1: float = element_parameters[2]
    beta1: float = element_parameters[3]
    gap: float = element_parameters[4]
    fint: float = element_parameters[5]
    r1: float = element_parameters[6]
    d: float = element_parameters[len(element_parameters) - 1]

    h = alpha / L
    h = h / (1 + d)
    k1 = k1 / (1 + d)
    T = np.zeros((6, 6, 6))
    phi1 = fint * h * gap * cos(beta1) ** (-1) * (1 + sin(beta1) ** 2)
    T[0, 0, 0] = -h * tan(beta1) ** 2 / 2
    T[0, 2, 2] = h * cos(beta1) ** (-1) ** 2 / 2
    T[1, 0, 0] = k1 * tan(beta1) + (h * cos(beta1) ** (-1)) ** 3 / (2 * r1)
    T[1, 0, 1] = h * tan(beta1) ** 2
    T[1, 0, 4] = -h * tan(beta1)
    T[1, 2, 2] = -k1 * tan(beta1) + 0.5 * h ** 2 * tan(beta1) + h ** 2 + tan(beta1) ** 3 - h * cos(
        beta1) ** (-1) ** 3 / (2 * r1)
    T[1, 2, 3] = -h * tan(beta1) ** 2
    T[2, 0, 2] = h * tan(beta1) ** 2
    T[3, 0, 2] = -2 * k1 * tan(beta1) - h * cos(beta1) ** (-1) ** 3 / r1
    T[3, 0, 3] = -h * tan(beta1) ** 2
    T[3, 1, 2] = -h * cos(beta1) ** (-1) ** 2
    T[3, 2, 4] = -h * phi1 * cos(beta1 - phi1) ** (-2) + h * tan(beta1)
    return T
