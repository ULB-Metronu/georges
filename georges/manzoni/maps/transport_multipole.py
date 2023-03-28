import numpy as np
from numba import njit
from numba.typed import List as nList
from numpy import cos, cosh, sin, sinh, sqrt


@njit(cache=True)
def compute_transport_multipole_matrix(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    k1: float = element_parameters[1]
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
def compute_transport_multipole_tensor(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    k1: float = element_parameters[1]
    k2: float = element_parameters[2]
    T = np.zeros((6, 6, 6))
    if k1 == 0:
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
    if k1 > 0:
        T[0, 0, 0] = k2 * (cos(L * sqrt(k1)) ** 2 + cos(L * sqrt(k1)) - 2) / (3 * k1)
        T[0, 0, 1] = 2 * k2 * (cos(L * sqrt(k1)) - 1) * sin(L * sqrt(k1)) / (3 * k1 ** (3 / 2))
        T[0, 0, 5] = L * sqrt(k1) * sin(L * sqrt(k1)) / 2
        T[0, 1, 1] = k2 * (sin(L * sqrt(k1)) ** 2 + 2 * cos(L * sqrt(k1)) - 2) / (3 * k1**2)
        T[0, 1, 5] = -L * cos(L * sqrt(k1)) / 2 + sin(L * sqrt(k1)) / (2 * sqrt(k1))
        T[0, 2, 2] = k2 * (-3 * cos(L * sqrt(k1)) + cosh(L * sqrt(k1)) ** 2 + 2) / (5 * k1)
        T[0, 2, 3] = -k2 * (2 * sin(L * sqrt(k1)) - sinh(2 * L * sqrt(k1))) / (5 * k1 ** (3 / 2))
        T[0, 3, 3] = k2 * (2 * cos(L * sqrt(k1)) + cosh(L * sqrt(k1)) ** 2 - 3) / (5 * k1**2)
        T[1, 0, 0] = -k2 * (sin(L * sqrt(k1)) + sin(2 * L * sqrt(k1))) / (3 * sqrt(k1))
        T[1, 0, 1] = 2 * k2 * (-cos(L * sqrt(k1)) + cos(2 * L * sqrt(k1))) / (3 * k1)
        T[1, 0, 5] = L * k1 * cos(L * sqrt(k1)) / 2 + sqrt(k1) * sin(L * sqrt(k1)) / 2
        T[1, 1, 1] = 2 * k2 * (cos(L * sqrt(k1)) - 1) * sin(L * sqrt(k1)) / (3 * k1 ** (3 / 2))
        T[1, 1, 5] = L * sqrt(k1) * sin(L * sqrt(k1)) / 2
        T[1, 2, 2] = k2 * (3 * sin(L * sqrt(k1)) + sinh(2 * L * sqrt(k1))) / (5 * sqrt(k1))
        T[1, 2, 3] = -2 * k2 * (cos(L * sqrt(k1)) - cosh(2 * L * sqrt(k1))) / (5 * k1)
        T[1, 3, 3] = -k2 * (2 * sin(L * sqrt(k1)) - sinh(2 * L * sqrt(k1))) / (5 * k1 ** (3 / 2))
        T[2, 0, 2] = (
            2
            * k2
            * (2 * sin(L * sqrt(k1)) * sinh(L * sqrt(k1)) - cos(L * sqrt(k1)) * cosh(L * sqrt(k1)) + cosh(L * sqrt(k1)))
            / (5 * k1)
        )
        T[2, 0, 3] = (
            -2
            * k2
            * (
                -2 * sin(L * sqrt(k1)) * cosh(L * sqrt(k1))
                + cos(L * sqrt(k1)) * sinh(L * sqrt(k1))
                + sinh(L * sqrt(k1))
            )
            / (5 * k1 ** (3 / 2))
        )
        T[2, 1, 2] = (
            -2
            * k2
            * (
                sin(L * sqrt(k1)) * cosh(L * sqrt(k1))
                + 2 * cos(L * sqrt(k1)) * sinh(L * sqrt(k1))
                - 3 * sinh(L * sqrt(k1))
            )
            / (5 * k1 ** (3 / 2))
        )
        T[2, 1, 3] = (
            -2
            * k2
            * (
                sin(L * sqrt(k1)) * sinh(L * sqrt(k1))
                + 2 * cos(L * sqrt(k1)) * cosh(L * sqrt(k1))
                - 2 * cosh(L * sqrt(k1))
            )
            / (5 * k1**2)
        )
        T[2, 2, 5] = -L * sqrt(k1) * sinh(L * sqrt(k1)) / 2
        T[2, 3, 5] = -L * cosh(L * sqrt(k1)) / 2 + sinh(L * sqrt(k1)) / (2 * sqrt(k1))
        T[3, 0, 2] = (
            2
            * k2
            * (3 * sin(L * sqrt(k1)) * cosh(L * sqrt(k1)) + cos(L * sqrt(k1)) * sinh(L * sqrt(k1)) + sinh(L * sqrt(k1)))
            / (5 * sqrt(k1))
        )
        T[3, 0, 3] = (
            2
            * k2
            * (3 * sin(L * sqrt(k1)) * sinh(L * sqrt(k1)) + cos(L * sqrt(k1)) * cosh(L * sqrt(k1)) - cosh(L * sqrt(k1)))
            / (5 * k1)
        )
        T[3, 1, 2] = (
            2
            * k2
            * (
                sin(L * sqrt(k1)) * sinh(L * sqrt(k1))
                - 3 * cos(L * sqrt(k1)) * cosh(L * sqrt(k1))
                + 3 * cosh(L * sqrt(k1))
            )
            / (5 * k1)
        )
        T[3, 1, 3] = (
            2
            * k2
            * (
                sin(L * sqrt(k1)) * cosh(L * sqrt(k1))
                - 3 * cos(L * sqrt(k1)) * sinh(L * sqrt(k1))
                + 2 * sinh(L * sqrt(k1))
            )
            / (5 * k1 ** (3 / 2))
        )
        T[3, 2, 5] = -L * k1 * cosh(L * sqrt(k1)) / 2 - sqrt(k1) * sinh(L * sqrt(k1)) / 2
        T[3, 3, 5] = -L * sqrt(k1) * sinh(L * sqrt(k1)) / 2
    if k1 < 0:
        T[0, 0, 0] = k2 * (cosh(L * sqrt(-k1)) ** 2 + cosh(L * sqrt(-k1)) - 2) / (3 * k1)
        T[0, 0, 1] = -2 * k2 * (cosh(L * sqrt(-k1)) - 1) * sinh(L * sqrt(-k1)) / (3 * (-k1) ** (3 / 2))
        T[0, 0, 5] = -L * sqrt(-k1) * sinh(L * sqrt(-k1)) / 2
        T[0, 1, 1] = -k2 * (sinh(L * sqrt(-k1)) ** 2 - 2 * cosh(L * sqrt(-k1)) + 2) / (3 * k1**2)
        T[0, 1, 5] = -L * cosh(L * sqrt(-k1)) / 2 + sinh(L * sqrt(-k1)) / (2 * sqrt(-k1))
        T[0, 2, 2] = -k2 * (sin(L * sqrt(-k1)) ** 2 + 3 * cosh(L * sqrt(-k1)) - 3) / (5 * k1)
        T[0, 2, 3] = -k2 * (sin(2 * L * sqrt(-k1)) - 2 * sinh(L * sqrt(-k1))) / (5 * (-k1) ** (3 / 2))
        T[0, 3, 3] = -k2 * (sin(L * sqrt(-k1)) ** 2 - 2 * cosh(L * sqrt(-k1)) + 2) / (5 * k1**2)
        T[1, 0, 0] = -k2 * (sinh(L * sqrt(-k1)) + sinh(2 * L * sqrt(-k1))) / (3 * sqrt(-k1))
        T[1, 0, 1] = -2 * k2 * (cosh(L * sqrt(-k1)) - cosh(2 * L * sqrt(-k1))) / (3 * k1)
        T[1, 0, 5] = L * k1 * cosh(L * sqrt(-k1)) / 2 - sqrt(-k1) * sinh(L * sqrt(-k1)) / 2
        T[1, 1, 1] = -2 * k2 * (cosh(L * sqrt(-k1)) - 1) * sinh(L * sqrt(-k1)) / (3 * (-k1) ** (3 / 2))
        T[1, 1, 5] = -L * sqrt(-k1) * sinh(L * sqrt(-k1)) / 2
        T[1, 2, 2] = k2 * (sin(2 * L * sqrt(-k1)) + 3 * sinh(L * sqrt(-k1))) / (5 * sqrt(-k1))
        T[1, 2, 3] = -2 * k2 * (2 * sin(L * sqrt(-k1)) ** 2 + cosh(L * sqrt(-k1)) - 1) / (5 * k1)
        T[1, 3, 3] = -k2 * (sin(2 * L * sqrt(-k1)) - 2 * sinh(L * sqrt(-k1))) / (5 * (-k1) ** (3 / 2))
        T[2, 0, 2] = (
            -2
            * k2
            * (
                2 * sin(L * sqrt(-k1)) * sinh(L * sqrt(-k1))
                + cos(L * sqrt(-k1)) * cosh(L * sqrt(-k1))
                - cos(L * sqrt(-k1))
            )
            / (5 * k1)
        )
        T[2, 0, 3] = (
            2
            * k2
            * (
                sin(L * sqrt(-k1)) * cosh(L * sqrt(-k1))
                + sin(L * sqrt(-k1))
                - 2 * cos(L * sqrt(-k1)) * sinh(L * sqrt(-k1))
            )
            / (5 * (-k1) ** (3 / 2))
        )
        T[2, 1, 2] = (
            2
            * k2
            * (
                2 * sin(L * sqrt(-k1)) * cosh(L * sqrt(-k1))
                - 3 * sin(L * sqrt(-k1))
                + cos(L * sqrt(-k1)) * sinh(L * sqrt(-k1))
            )
            / (5 * (-k1) ** (3 / 2))
        )
        T[2, 1, 3] = (
            2
            * k2
            * (
                sin(L * sqrt(-k1)) * sinh(L * sqrt(-k1))
                - 2 * cos(L * sqrt(-k1)) * cosh(L * sqrt(-k1))
                + 2 * cos(L * sqrt(-k1))
            )
            / (5 * k1**2)
        )
        T[2, 2, 5] = L * sqrt(-k1) * sin(L * sqrt(-k1)) / 2
        T[2, 3, 5] = -L * cos(L * sqrt(-k1)) / 2 + sin(L * sqrt(-k1)) / (2 * sqrt(-k1))
        T[3, 0, 2] = (
            2
            * k2
            * (
                sin(L * sqrt(-k1)) * cosh(L * sqrt(-k1))
                + sin(L * sqrt(-k1))
                + 3 * cos(L * sqrt(-k1)) * sinh(L * sqrt(-k1))
            )
            / (5 * sqrt(-k1))
        )
        T[3, 0, 3] = (
            -2
            * k2
            * (
                3 * sin(L * sqrt(-k1)) * sinh(L * sqrt(-k1))
                - cos(L * sqrt(-k1)) * cosh(L * sqrt(-k1))
                + cos(L * sqrt(-k1))
            )
            / (5 * k1)
        )
        T[3, 1, 2] = (
            -2
            * k2
            * (
                sin(L * sqrt(-k1)) * sinh(L * sqrt(-k1))
                + 3 * cos(L * sqrt(-k1)) * cosh(L * sqrt(-k1))
                - 3 * cos(L * sqrt(-k1))
            )
            / (5 * k1)
        )
        T[3, 1, 3] = (
            -2
            * k2
            * (
                -3 * sin(L * sqrt(-k1)) * cosh(L * sqrt(-k1))
                + 2 * sin(L * sqrt(-k1))
                + cos(L * sqrt(-k1)) * sinh(L * sqrt(-k1))
            )
            / (5 * (-k1) ** (3 / 2))
        )
        T[3, 2, 5] = -L * k1 * cos(L * sqrt(-k1)) / 2 + sqrt(-k1) * sin(L * sqrt(-k1)) / 2
        T[3, 3, 5] = L * sqrt(-k1) * sin(L * sqrt(-k1)) / 2
    return T
