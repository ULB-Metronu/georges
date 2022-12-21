import numpy as np
from numba import njit
from numba.typed import List as nList
from numpy import cos, cosh, sin, sinh, sqrt


@njit(cache=True)
def compute_transport_combined_dipole_ex_matrix(
    element_parameters: nList,
) -> np.ndarray:

    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    h = alpha / L
    k1: float = element_parameters[2]
    d: float = element_parameters[len(element_parameters) - 1]
    h = h / (1 + d)
    k1 = k1 / (1 + d)

    R = np.zeros((6, 6))
    R[4, 4] = 1
    R[5, 5] = 1
    if h == 0:
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
    if h != 0:
        if k1 == 0:
            R[0, 0] = cos(L * h)
            R[0, 1] = sin(L * h) / h
            R[0, 5] = (1 - cos(L * h)) / h
            R[1, 0] = -h * sin(L * h)
            R[1, 1] = cos(L * h)
            R[1, 5] = sin(L * h)
            R[2, 2] = 1
            R[2, 3] = L
            R[3, 3] = 1
        if k1 != 0:
            if h**2 + k1 == 0:
                R[0, 0] = 1
                R[0, 1] = L
                R[0, 5] = L**2 * h / 2
                R[1, 1] = 1
                R[1, 5] = L * h
                R[2, 2] = cos(L * sqrt(-k1))
                R[2, 3] = sin(L * sqrt(-k1)) / sqrt(-k1)
                R[3, 2] = -sqrt(-k1) * sin(L * sqrt(-k1))
                R[3, 3] = cos(L * sqrt(-k1))
            if h**2 + k1 > 0:
                if k1 > 0:
                    R[0, 0] = cos(L * sqrt(h**2 + k1))
                    R[0, 1] = sin(L * sqrt(h**2 + k1)) / sqrt(h**2 + k1)
                    R[0, 5] = h * (1 - cos(L * sqrt(h**2 + k1))) / (h**2 + k1)
                    R[1, 0] = -sqrt(h**2 + k1) * sin(L * sqrt(h**2 + k1))
                    R[1, 1] = cos(L * sqrt(h**2 + k1))
                    R[1, 5] = h * sin(L * sqrt(h**2 + k1)) / sqrt(h**2 + k1)
                    R[2, 2] = cosh(L * sqrt(k1))
                    R[2, 3] = sinh(L * sqrt(k1)) / sqrt(k1)
                    R[3, 2] = sqrt(k1) * sinh(L * sqrt(k1))
                    R[3, 3] = cosh(L * sqrt(k1))
                if k1 < 0:
                    R[0, 0] = cos(L * sqrt(h**2 + k1))
                    R[0, 1] = sin(L * sqrt(h**2 + k1)) / sqrt(h**2 + k1)
                    R[0, 5] = h * (1 - cos(L * sqrt(h**2 + k1))) / (h**2 + k1)
                    R[1, 0] = -sqrt(h**2 + k1) * sin(L * sqrt(h**2 + k1))
                    R[1, 1] = cos(L * sqrt(h**2 + k1))
                    R[1, 5] = h * sin(L * sqrt(h**2 + k1)) / sqrt(h**2 + k1)
                    R[2, 2] = cos(L * sqrt(-k1))
                    R[2, 3] = sin(L * sqrt(-k1)) / sqrt(-k1)
                    R[3, 2] = -sqrt(-k1) * sin(L * sqrt(-k1))
                    R[3, 3] = cos(L * sqrt(-k1))
            if h**2 + k1 < 0:
                R[0, 0] = cosh(L * sqrt(-(h**2) - k1))
                R[0, 1] = sinh(L * sqrt(-(h**2) - k1)) / sqrt(-(h**2) - k1)
                R[0, 5] = -(h * cosh(L * sqrt(-(h**2) - k1)) - h) / (h**2 + k1)
                R[1, 0] = sqrt(-(h**2) - k1) * sinh(L * sqrt(-(h**2) - k1))
                R[1, 1] = cosh(L * sqrt(-(h**2) - k1))
                R[1, 5] = h * sinh(L * sqrt(-(h**2) - k1)) / sqrt(-(h**2) - k1)
                R[2, 2] = cos(L * sqrt(-k1))
                R[2, 3] = sin(L * sqrt(-k1)) / sqrt(-k1)
                R[3, 2] = -sqrt(-k1) * sin(L * sqrt(-k1))
                R[3, 3] = cos(L * sqrt(-k1))
    return R


@njit(cache=True)
def compute_transport_combined_dipole_ex_tensor(
    element_parameters: nList,
) -> np.ndarray:

    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    h = alpha / L
    k1: float = element_parameters[2]
    d: float = element_parameters[len(element_parameters) - 1]
    h = h / (1 + d)
    k1 = k1 / (1 + d)

    T = np.zeros((6, 6, 6))
    if h == 0:
        if k1 > 0:
            T[0, 0, 5] = L**2 * k1 / 2
            T[0, 1, 5] = L**3 * k1 / 6
            T[1, 0, 5] = L * k1
            T[1, 1, 5] = L**2 * k1 / 2
            T[2, 2, 5] = -L * sqrt(k1) * sinh(L * sqrt(k1)) / 2
            T[2, 3, 5] = -L * cosh(L * sqrt(k1)) / 2 + sinh(
                L * sqrt(k1),
            ) / (2 * sqrt(k1))
            T[3, 2, 5] = (
                -L * k1 * cosh(L * sqrt(k1)) / 2
                - sqrt(k1)
                * sinh(
                    L * sqrt(k1),
                )
                / 2
            )
            T[3, 3, 5] = -L * sqrt(k1) * sinh(L * sqrt(k1)) / 2
        if k1 < 0:
            T[0, 0, 5] = L**2 * k1 / 2
            T[0, 1, 5] = L**3 * k1 / 6
            T[1, 0, 5] = L * k1
            T[1, 1, 5] = L**2 * k1 / 2
            T[2, 2, 5] = L * sqrt(-k1) * sin(L * sqrt(-k1)) / 2
            T[2, 3, 5] = -L * cos(L * sqrt(-k1)) / 2 + sin(
                L * sqrt(-k1),
            ) / (2 * sqrt(-k1))
            T[3, 2, 5] = (
                -L * k1 * cos(L * sqrt(-k1)) / 2
                + sqrt(-k1)
                * sin(
                    L * sqrt(-k1),
                )
                / 2
            )
            T[3, 3, 5] = L * sqrt(-k1) * sin(L * sqrt(-k1)) / 2
    if h != 0:
        if k1 == 0:
            T[0, 0, 0] = -h * sin(L * h) ** 2 / 2
            T[0, 0, 1] = (cos(L * h) - 1) * sin(L * h)
            T[0, 0, 5] = 1 - cos(L * h) ** 2
            T[0, 1, 1] = (1 - cos(L * h)) * cos(L * h) / (2 * h)
            T[0, 1, 5] = (-12 * L * h * cos(L * h) + 10 * sin(L * h) + sin(2 * L * h)) / (6 * h)
            T[0, 3, 3] = (cos(L * h) - 1) / (2 * h)
            T[0, 5, 5] = (
                -12 * L * h * sin(L * h)
                + 6 * sin(L * h) ** 6
                - 18 * sin(L * h) ** 4
                + 19 * sin(L * h) ** 2
                + 6 * cos(L * h) ** 6
                - 16 * cos(L * h)
                + 10
            ) / (6 * h)
            T[1, 0, 0] = -(h**2) * sin(2 * L * h) / 2
            T[1, 0, 1] = h * (-cos(L * h) + cos(2 * L * h))
            T[1, 0, 5] = h * sin(2 * L * h)
            T[1, 1, 1] = -sin(L * h) / 2 + sin(2 * L * h) / 2
            T[1, 1, 5] = (
                2 * L * h * sin(L * h)
                + 2 * cos(L * h) ** 2 / 3
                - cos(
                    L * h,
                )
                / 3
                - 1 / 3
            )
            T[1, 3, 3] = -sin(L * h) / 2
            T[1, 5, 5] = (
                -2 * L * h * cos(L * h)
                + sin(L * h)
                * cos(
                    L * h,
                )
                / 3
                + 2 * sin(L * h) / 3
            )
            T[2, 0, 3] = -L * h + sin(L * h)
            T[2, 1, 3] = (1 - cos(L * h)) / h
            T[2, 3, 5] = L - sin(L * h) / h
            T[3, 0, 3] = h * (cos(L * h) - 1)
            T[3, 1, 3] = sin(L * h)
            T[3, 3, 5] = 1 - cos(L * h)
    return T
