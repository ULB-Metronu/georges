from numba import njit
import numpy as np
from numba.typed import List as nList


@njit
def compute_transport_fringe_out_Ex_matrix(element_parameters: nList) -> np.ndarray:
    h: float = element_parameters[0]
    beta2: float = element_parameters[3]
    fintx: float = element_parameters[6]
    gap: float = element_parameters[5]
    d: float = element_parameters[8]
    R = np.zeros((6, 6))
    h = h / (1 + d)
    phi2 = fintx * h * gap * np.cos(beta2) ** (-1) * (1 + np.sin(beta2) ** 2)
    R[0, 0] = 1
    R[1, 0] = h * np.tan(beta2)
    R[1, 1] = 1
    R[2, 2] = 1
    R[3, 2] = -h * np.tan(beta2 - phi2)
    R[3, 3] = 1
    R[4, 4] = 1
    R[5, 5] = 1
    return R


def compute_transport_fringe_out_Ex_tensor(element_parameters: nList) -> np.ndarray:
    h: float = element_parameters[0]
    k1: float = element_parameters[1]
    beta2: float = element_parameters[3]
    r2: float = element_parameters[7]
    fintx: float = element_parameters[6]
    gap: float = element_parameters[5]
    d: float = element_parameters[8]
    T = np.zeros((6, 6, 6))
    h = h / (1 + d)
    k1 = k1 / (1 + d)
    phi2 = fintx * h * gap * np.cos(beta2) ** (-1) * (1 + np.sin(beta2) ** 2)
    T[0, 0, 0] = h * np.tan(beta2) ** 2 / 2
    T[0, 2, 2] = -h * np.cos(beta2) ** (-1) ** 2 / 2
    T[1, 0, 0] = k1 * np.tan(beta2) - 0.5 * h ** 2 * np.tan(beta2) ** 3 + h * np.cos(beta2) ** (-1) ** 3 / (2 * r2)
    T[1, 0, 1] = -h * np.tan(beta2) ** 2
    T[1, 0, 4] = -h * np.tan(beta2)
    T[1, 2, 2] = -k1 * np.tan(beta2) - 0.5 * h ** 2 * np.tan(beta2) ** 3 - h * np.cos(beta2) ** (-1) ** 3 / (2 * r2)
    T[1, 2, 3] = h * np.tan(beta2) ** 2
    T[2, 0, 2] = -h * np.tan(beta2) ** 2
    T[3, 0, 2] = -2 * k1 * np.tan(beta2) + h ** 2 * np.tan(beta2) * np.cos(beta2) ** (-1) ** 2 - \
                 h * np.cos(beta2) ** (-1) ** 3 / r2
    T[3, 0, 3] = h * np.tan(beta2) ** 2
    T[3, 1, 2] = h * np.cos(beta2) ** (-1) ** 2
    T[3, 2, 4] = -h * phi2 * np.cos(beta2 - phi2) ** (-2) + h * np.tan(beta2)
    return T
