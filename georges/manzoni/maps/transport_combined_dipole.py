from numba import njit
import numpy as np
from numba.typed import List as nList


@njit
def compute_transport_combined_dipole_matrix(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    h = alpha / L
    k1: float = element_parameters[2]
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
            R[0, 0] = np.cos(L * np.sqrt(k1))
            R[0, 1] = np.sin(L * np.sqrt(k1)) / np.sqrt(k1)
            R[1, 0] = -np.sqrt(k1) * np.sin(L * np.sqrt(k1))
            R[1, 1] = np.cos(L * np.sqrt(k1))
            R[2, 2] = np.cosh(L * np.sqrt(k1))
            R[2, 3] = np.sinh(L * np.sqrt(k1)) / np.sqrt(k1)
            R[3, 2] = np.sqrt(k1) * np.sinh(L * np.sqrt(k1))
            R[3, 3] = np.cosh(L * np.sqrt(k1))
        if k1 < 0:
            R[0, 0] = np.cosh(L * np.sqrt(-k1))
            R[0, 1] = np.sinh(L * np.sqrt(-k1)) / np.sqrt(-k1)
            R[1, 0] = np.sqrt(-k1) * np.sinh(L * np.sqrt(-k1))
            R[1, 1] = np.cosh(L * np.sqrt(-k1))
            R[2, 2] = np.cos(L * np.sqrt(-k1))
            R[2, 3] = np.sin(L * np.sqrt(-k1)) / np.sqrt(-k1)
            R[3, 2] = -np.sqrt(-k1) * np.sin(L * np.sqrt(-k1))
            R[3, 3] = np.cos(L * np.sqrt(-k1))
    if h != 0:
        if k1 == 0:
            R[0, 0] = np.cos(L * h)
            R[0, 1] = np.sin(L * h) / h
            R[0, 4] = (1 - np.cos(L * h)) / h
            R[1, 0] = -np.sin(L * h) * h
            R[1, 1] = np.cos(L * h)
            R[2, 2] = 1
            R[2, 3] = L
            R[3, 3] = 1
        if k1 != 0:
            if h ** 2 + k1 == 0:
                R[0, 0] = 1
                R[0, 1] = L
                R[0, 4] = L ** 2 * h / 2
                R[1, 1] = 1
                R[2, 2] = np.cos(L * np.sqrt(-k1))
                R[2, 3] = np.sin(L * np.sqrt(-k1)) / np.sqrt(-k1)
                R[3, 2] = -np.sqrt(-k1) * np.sin(L * np.sqrt(-k1))
                R[3, 3] = np.cos(L * np.sqrt(-k1))
            if h ** 2 + k1 > 0:
                if k1 > 0:
                    R[0, 0] = np.cos(L * np.sqrt(k1 + h ** 2))
                    R[0, 1] = np.sin(L * np.sqrt(k1 + h ** 2)) / np.sqrt(k1 + h ** 2)
                    R[0, 4] = h * (1 - np.cos(L * np.sqrt(k1 + h ** 2))) / (k1 + h ** 2)
                    R[1, 0] = -np.sqrt(k1 + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2))
                    R[1, 1] = np.cos(L * np.sqrt(k1 + h ** 2))
                    R[2, 2] = np.cosh(np.sqrt(k1) * L)
                    R[2, 3] = np.sinh(np.sqrt(k1) * L) / np.sqrt(k1)
                    R[3, 2] = np.sqrt(k1) * np.sinh(np.sqrt(k1) * L)
                    R[3, 3] = np.cosh(np.sqrt(k1) * L)
                if k1 < 0:
                    R[0, 0] = np.cos(L * np.sqrt(k1 + h ** 2))
                    R[0, 1] = np.sin(L * np.sqrt(k1 + h ** 2)) / np.sqrt(k1 + h ** 2)
                    R[0, 4] = h * (1 - np.cos(L * np.sqrt(k1 + h ** 2))) / (k1 + h ** 2)
                    R[1, 0] = -np.sqrt(k1 + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2))
                    R[1, 1] = np.cos(L * np.sqrt(k1 + h ** 2))
                    R[2, 2] = np.cos(L * np.sqrt(-k1))
                    R[2, 3] = np.sin(L * np.sqrt(-k1)) / np.sqrt(-k1)
                    R[3, 2] = -np.sqrt(-k1) * np.sin(L * np.sqrt(-k1))
                    R[3, 3] = np.cos(L * np.sqrt(-k1))
            if h ** 2 + k1 < 0:
                R[0, 0] = np.cosh(L * np.sqrt(-k1 - h ** 2))
                R[0, 1] = np.sinh(L * np.sqrt(-k1 - h ** 2)) / np.sqrt(-k1 - h ** 2)
                R[0, 4] = h * (1 - np.cosh(L * np.sqrt(-k1 - h ** 2))) / (k1 + h ** 2)
                R[1, 0] = np.sqrt(-k1 - h ** 2) * np.sinh(L * np.sqrt(-k1 - h ** 2))
                R[1, 1] = np.cosh(L * np.sqrt(-k1 - h ** 2))
                R[2, 2] = np.cos(L * np.sqrt(-k1))
                R[2, 3] = np.sin(L * np.sqrt(-k1)) / np.sqrt(-k1)
                R[3, 2] = -np.sqrt(-k1) * np.sin(L * np.sqrt(-k1))
                R[3, 3] = np.cos(L * np.sqrt(-k1))
    return R


@njit
def compute_transport_combined_dipole_tensor(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    h = alpha / L
    k1: float = element_parameters[2]
    T = np.zeros((6, 6, 6))
    if h == 0:
        if k1 > 0:
            T[0, 0, 4] = L * np.sqrt(k1) * np.sin(L * np.sqrt(k1)) / 2
            T[0, 1, 4] = (-L * k1 * np.cos(L * np.sqrt(k1)) + np.sqrt(k1) * np.sin(L * np.sqrt(k1))) / (2 * k1)
            T[1, 0, 4] = L * k1 * np.cos(L * np.sqrt(k1)) / 2 + np.sqrt(k1) * np.sin(L * np.sqrt(k1)) / 2
            T[1, 1, 4] = L * np.sqrt(k1) * np.sin(L * np.sqrt(k1)) / 2
            T[2, 2, 4] = -L * np.sqrt(k1) * np.sinh(L * np.sqrt(k1)) / 2
            T[2, 3, 4] = -L * np.cosh(L * np.sqrt(k1)) / 2 + np.sinh(L * np.sqrt(k1)) / (2 * np.sqrt(k1))
            T[3, 2, 4] = -L * k1 * np.cosh(L * np.sqrt(k1)) / 2 - np.sqrt(k1) * np.sinh(L * np.sqrt(k1)) / 2
            T[3, 3, 4] = -L * np.sqrt(k1) * np.sinh(L * np.sqrt(k1)) / 2
        if k1 < 0:
            T[0, 0, 4] = L * k1 * np.sinh(L * np.sqrt(-k1)) / (2 * np.sqrt(-k1))
            T[0, 1, 4] = L * np.cosh(L * np.sqrt(-k1)) / 2 + k1 * np.sinh(L * np.sqrt(-k1)) / (2 * (-k1) ** (3 / 2))
            T[1, 0, 4] = L * k1 * np.cosh(L * np.sqrt(-k1)) / 2 + k1 * np.sinh(L * np.sqrt(-k1)) / (2 * np.sqrt(-k1))
            T[1, 1, 4] = -L * k1 * np.sinh(L * np.sqrt(-k1)) / (2 * np.sqrt(-k1))
            T[2, 2, 4] = L * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) / 2
            T[2, 3, 4] = L * np.cos(L * np.sqrt(-k1)) / 2 - np.sin(L * np.sqrt(-k1)) / (2 * np.sqrt(-k1))
            T[3, 2, 4] = -L * k1 * np.cos(L * np.sqrt(-k1)) / 2 + np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) / 2
            T[3, 3, 4] = -L * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) / 2
    if h != 0:
        if k1 == 0:
            T[0, 0, 0] = (3 * h ** 2 * np.cos(L * h) ** 2 - 3 * h ** 2) / (6 * h)
            T[0, 0, 1] = (np.cos(L * h) - 1) * np.sin(L * h) * h / h
            T[0, 0, 4] = (h ** 8 * (np.sin(L * h) ** 2 + 2 * np.cos(L * h) - 2) / 3 + h ** 8 * (
                    -3 * np.sin(L * h) ** 4 * np.cos(L * h) + 4 * np.sin(L * h) ** 2 + 3 * np.cos(L * h) ** 5 - 6 * np.cos(
                L * h) ** 3 - np.cos(L * h) + 4) / 6) / h ** 8
            T[0, 1, 1] = (h ** 2 * (np.sin(L * h) ** 2 - np.cos(L * h) + 1) + 2 * h ** 2 * (
                    np.sin(L * h) ** 2 + 2 * np.cos(L * h) - 2)) / (6 * h ** 3)
            T[0, 1, 4] = (-2 * L * h ** 12 * np.cos(L * h) + h ** 10 * (4 - np.cos(L * h)) * np.sin(L * h) * h / 3 + h ** 10 * (
                    np.sin(L * h) + np.sin(2 * L * h)) * h / 3) / h ** 12
            T[0, 3, 3] = (np.cos(L * h) - 1) / (2 * h)
            T[0, 4, 4] = (-12 * L * h ** 18 * np.sin(L * h) * h + 6 * h ** 18 * (np.cos(L * h) - 1) + h ** 18 * (
                    np.cos(L * h) ** 2 - 14 * np.cos(L * h) + 13) + h ** 18 * (
                                  6 * np.sin(L * h) ** 6 - 3 * np.sin(L * h) ** 4 * np.cos(L * h) - 18 * np.sin(
                              L * h) ** 4 + 20 * np.sin(L * h) ** 2 + 6 * np.cos(L * h) ** 6 + 3 * np.cos(
                              L * h) ** 5 - 6 * np.cos(L * h) ** 3 - 5 * np.cos(L * h) + 2)) / (6 * h ** 19)
            T[1, 0, 0] = -h ** 3 * np.sin(L * h) * np.cos(L * h) / h
            T[1, 0, 1] = h * (-np.cos(L * h) + np.cos(2 * L * h))
            T[1, 0, 4] = np.sin(2 * L * h) * h
            T[1, 1, 1] = (4 * h ** 2 * (np.cos(L * h) - 1) + h ** 2 * (2 * np.cos(L * h) + 1)) * np.sin(L * h) * h / (6 * h ** 3)
            T[1, 1, 4] = (6 * L * h ** 12 * np.sin(L * h) * h + 2 * h ** 12 * np.cos(L * h) ** 2 - h ** 12 * np.cos(
                L * h) - h ** 12) / (3 * h ** 12)
            T[1, 3, 3] = -h * np.sin(L * h) / (2 * h)
            T[1, 4, 4] = -2 * L * h * np.cos(L * h) - h * np.sin(L * h) / h + 5 * np.sin(L * h) * h / (3 * h) + np.sin(
                2 * L * h) * h / (6 * h)
            T[2, 0, 3] = -L * h + h * np.sin(L * h) / h
            T[2, 1, 3] = (1 - np.cos(L * h)) / h
            T[2, 3, 4] = L - np.sin(L * h) * h / h ** 2
            T[3, 0, 3] = h * (np.cos(L * h) - 1)
            T[3, 1, 3] = h * np.sin(L * h) / h
            T[3, 3, 4] = 1 - np.cos(L * h)
        if k1 != 0:
            if h ** 2 + k1 == 0:
                T[0, 0, 0] = -L ** 2 * h * (2 * k1 + h ** 2) / 2
                T[0, 0, 1] = -L ** 3 * h * (2 * k1 + h ** 2) / 3
                T[0, 0, 4] = L ** 2 * (-k1 * L ** 2 * h ** 2 / 6 + k1 / 2 - L ** 2 * h ** 4 / 12 + h ** 2)
                T[0, 1, 1] = L ** 2 * h * (-2 * k1 * L ** 2 - L ** 2 * h ** 2 + 3) / 12
                T[0, 1, 4] = L ** 3 * (k1 * L ** 2 * h ** 2 / 10 + k1 / 6 + L ** 2 * h ** 4 / 20 + h ** 2 / 2)
                T[0, 2, 2] = h * np.sin(L * np.sqrt(-k1)) ** 2 / 4
                T[0, 3, 3] = -L ** 2 * h / 4
                T[0, 4, 4] = L ** 2 * h * (
                        2 * k1 * L ** 4 * h ** 2 + 5 * k1 * L ** 2 + L ** 4 * h ** 4 + 15 * L ** 2 * h ** 2 - 60) / 120
                T[1, 0, 0] = -L * h * (2 * k1 + h ** 2)
                T[1, 0, 1] = -L ** 2 * h * (2 * k1 + h ** 2)
                T[1, 0, 4] = L * (-2 * k1 * L ** 2 * h ** 2 + 3 * k1 - L ** 2 * h ** 4 + 6 * h ** 2) / 3
                T[1, 1, 1] = L * h * (-4 * k1 * L ** 2 - 2 * L ** 2 * h ** 2 + 3) / 6
                T[1, 1, 4] = L ** 2 * (2 * k1 * L ** 2 * h ** 2 + 2 * k1 + L ** 2 * h ** 4 + 6 * h ** 2) / 4
                T[1, 2, 2] = h * np.sqrt(-k1) * np.sin(2 * L * np.sqrt(-k1)) / 4
                T[1, 3, 3] = -L * h / 2
                T[1, 4, 4] = L * h * (
                        k1 * L ** 4 * h ** 2 / 10 + k1 * L ** 2 / 6 + L ** 4 * h ** 4 / 20 + L ** 2 * h ** 2 / 2 - 1)
                T[2, 0, 2] = -L * h * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1))
                T[2, 0, 3] = -L * h * np.cos(L * np.sqrt(-k1)) + h * np.sin(L * np.sqrt(-k1)) / np.sqrt(-k1)
                T[2, 1, 2] = -L ** 2 * h * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) / 2
                T[2, 1, 3] = L ** 2 * h * np.cos(L * np.sqrt(-k1)) / 2
                T[2, 2, 4] = -L * np.sqrt(-k1) * (L ** 2 * h ** 2 - 3) * np.sin(L * np.sqrt(-k1)) / 6
                T[2, 3, 4] = -L ** 3 * h ** 2 * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) / 6 - L ** 2 * h ** 2 * np.cos(
                    L * np.sqrt(-k1)) / 4 + L ** 2 * h ** 2 * np.sin(L * np.sqrt(-k1)) / (4 * np.sqrt(-k1)) + L * h ** 2 * np.sin(
                    L * np.sqrt(-k1)) / (4 * np.sqrt(-k1)) - L * np.cos(L * np.sqrt(-k1)) / 2 - h ** 2 * np.sin(L * np.sqrt(-k1)) / (
                                     4 * (-k1) ** (3 / 2)) + np.sin(L * np.sqrt(-k1)) / (
                                     2 * np.sqrt(-k1)) - L * h ** 2 * np.cos(L * np.sqrt(-k1)) / (4 * k1)
                T[3, 0, 2] = h * (k1 * L * np.cos(L * np.sqrt(-k1)) - np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)))
                T[3, 0, 3] = L * h * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1))
                T[3, 1, 2] = L * h * (k1 * L * np.cos(L * np.sqrt(-k1)) - 2 * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1))) / 2
                T[3, 1, 3] = L * h * (-L * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) + 2 * np.cos(L * np.sqrt(-k1))) / 2
                T[3, 2, 4] = k1 * L * (L ** 2 * h ** 2 - 3) * np.cos(L * np.sqrt(-k1)) / 6 - L ** 2 * h ** 2 * np.sqrt(
                    -k1) * np.sin(L * np.sqrt(-k1)) / 3 - np.sqrt(-k1) * (L ** 2 * h ** 2 - 3) * np.sin(L * np.sqrt(-k1)) / 6
                T[3, 3, 4] = (3 * k1 * L * h ** 2 * np.sin(L * np.sqrt(-k1)) - 3 * k1 * h ** 2 * (2 * L + 1) * np.sin(
                    L * np.sqrt(-k1)) - L * (-k1) ** (3 / 2) * (
                                      -2 * k1 * L ** 2 * h ** 2 * np.cos(L * np.sqrt(-k1)) + 3 * L * h ** 2 * np.sqrt(
                                  -k1) * np.sin(L * np.sqrt(-k1)) - 3 * L * h ** 2 * np.cos(
                                  L * np.sqrt(-k1)) + 3 * h ** 2 * np.cos(L * np.sqrt(-k1)) - 6 * np.sqrt(-k1) * np.sin(
                                  L * np.sqrt(-k1)))) / (12 * (-k1) ** (3 / 2))
            if h ** 2 + k1 > 0:
                if k1 > 0:
                    T[0, 0, 0] = h * (4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) ** 2 + 4 * k1 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 7 * k1 + 2 * h ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) ** 2 + 2 * h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) - 3 * h ** 2 + (
                                              k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2)) ** 2 - 2 * (
                                              k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2))) / (6 * (k1 + h ** 2))
                    T[0, 0, 1] = h * (4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) - 5 * k1 + 2 * h ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 3 * h ** 2 + (k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2))) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) / (3 * (k1 + h ** 2) ** (3 / 2))
                    T[0, 0, 4] = (-6 * L * h ** 2 * (k1 + h ** 2) ** (7 / 2) * (2 * k1 + h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + 3 * L * (k1 + h ** 2) ** (9 / 2) * (k1 + 2 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + 2 * h ** 2 * (k1 + h ** 2) ** 4 * (
                                          np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 + 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) - 2) + h ** 2 * (k1 + h ** 2) ** 3 * (
                                          -6 * k1 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) + 8 * k1 * np.sin(
                                      L * np.sqrt(k1 + h ** 2)) ** 2 + 6 * k1 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) ** 5 - 12 * k1 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) ** 3 - 2 * k1 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) + 8 * k1 - 3 * h ** 2 * np.sin(
                                      L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(L * np.sqrt(k1 + h ** 2)) + 4 * h ** 2 * np.sin(
                                      L * np.sqrt(k1 + h ** 2)) ** 2 + 3 * h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) ** 5 - 6 * h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) ** 3 - h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) + 4 * h ** 2)) / (6 * (k1 + h ** 2) ** 5)
                    T[0, 1, 1] = h * (4 * k1 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 + 8 * k1 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 7 * k1 + 2 * h ** 2 * np.sin(
                        L * np.sqrt(k1 + h ** 2)) ** 2 + 4 * h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) - 3 * h ** 2 + (
                                              k1 + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 - (k1 + h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2))) / (6 * (k1 + h ** 2) ** 2)
                    T[0, 1, 4] = (-L * h ** 2 * (k1 + h ** 2) ** 5 * (2 * k1 + h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - L * (k1 + h ** 2) ** 6 * (k1 + 2 * h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) / 2 + h ** 2 * (k1 + h ** 2) ** (9 / 2) * (
                                          4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) + 2 * k1 + 2 * h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) / 3 + (
                                          k1 + h ** 2) ** (11 / 2) * (
                                          3 * k1 - 2 * h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) + 8 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) / 6) / (k1 + h ** 2) ** 7
                    T[0, 2, 2] = k1 * h * (np.cos(L * np.sqrt(k1 + h ** 2)) - np.cosh(2 * np.sqrt(k1) * L)) / (
                            2 * (5 * k1 + h ** 2))
                    T[0, 3, 3] = h * (np.cos(L * np.sqrt(k1 + h ** 2)) - 1) / (2 * (k1 + h ** 2))
                    T[0, 4, 4] = h * (-6 * L * h ** 2 * (k1 + h ** 2) ** (15 / 2) * (2 * k1 + h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) - 3 * L * (k1 + h ** 2) ** (17 / 2) * (k1 + 2 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + h ** 2 * (k1 + h ** 2) ** 7 * (
                                              12 * k1 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 6 - 6 * k1 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(L * np.sqrt(k1 + h ** 2)) - 36 * k1 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 4 + 40 * k1 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 2 + 12 * k1 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 6 + 6 * k1 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 5 - 12 * k1 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 3 - 10 * k1 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) + 4 * k1 + 6 * h ** 2 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 6 - 3 * h ** 2 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) - 18 * h ** 2 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 4 + 20 * h ** 2 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 2 + 6 * h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 6 + 3 * h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 5 - 6 * h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 3 - 5 * h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) + 2 * h ** 2) + 6 * (k1 + h ** 2) ** 9 * (
                                              np.cos(L * np.sqrt(k1 + h ** 2)) - 1) + (k1 + h ** 2) ** 8 * (
                                              -6 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) + 6 * k1 + h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 2 - 14 * h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) + 13 * h ** 2)) / (6 * (k1 + h ** 2) ** 10)
                    T[1, 0, 0] = -h * (
                            4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) + k1 + 2 * h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) + (
                            k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2))) * np.sin(L * np.sqrt(k1 + h ** 2)) / (
                                         3 * np.sqrt(k1 + h ** 2))
                    T[1, 0, 1] = h * (-(5 * k1 + 3 * h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 + (
                            4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) - 5 * k1 + 2 * h ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 3 * h ** 2 + (k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2))) * np.cos(
                        L * np.sqrt(k1 + h ** 2))) / (3 * (k1 + h ** 2))
                    T[1, 0, 4] = (-L * h ** 2 * (k1 + h ** 2) ** 4 * (2 * k1 + h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + L * (k1 + h ** 2) ** 5 * (k1 + 2 * h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) / 2 + 2 * h ** 2 * (k1 + h ** 2) ** (9 / 2) * (
                                          np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) / 3 - h ** 2 * (k1 + h ** 2) ** (7 / 2) * (2 * k1 + h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + 2 * h ** 2 * (k1 + h ** 2) ** (7 / 2) * (
                                          4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) + 2 * k1 + 2 * h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) / 3 + (
                                          k1 + h ** 2) ** (9 / 2) * (k1 + 2 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) / 2) / (k1 + h ** 2) ** 5
                    T[1, 1, 1] = h * (8 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) - 7 * k1 + 4 * h ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 3 * h ** 2 + 2 * (k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2))) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) / (6 * (k1 + h ** 2) ** (3 / 2))
                    T[1, 1, 4] = 2 * k1 * L * h ** 2 * np.sin(L * np.sqrt(k1 + h ** 2)) / (k1 + h ** 2) ** (
                            3 / 2) + k1 * L * np.sin(L * np.sqrt(k1 + h ** 2)) / (
                                         2 * np.sqrt(k1 + h ** 2)) - 8 * k1 * h ** 2 * np.sin(
                        L * np.sqrt(k1 + h ** 2)) ** 2 / (3 * (k1 + h ** 2) ** 2) - 4 * k1 * h ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) / (3 * (k1 + h ** 2) ** 2) + 4 * k1 * h ** 2 / (
                                         3 * (k1 + h ** 2) ** 2) + L * h ** 4 * np.sin(L * np.sqrt(k1 + h ** 2)) / (
                                         k1 + h ** 2) ** (3 / 2) + L * h ** 2 * np.sin(L * np.sqrt(k1 + h ** 2)) / np.sqrt(
                        k1 + h ** 2) - 4 * h ** 4 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 / (
                                         3 * (k1 + h ** 2) ** 2) - 2 * h ** 4 * np.cos(L * np.sqrt(k1 + h ** 2)) / (
                                         3 * (k1 + h ** 2) ** 2) + 2 * h ** 4 / (
                                         3 * (k1 + h ** 2) ** 2) + 2 * h ** 2 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 / (
                                         3 * (k1 + h ** 2)) + h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) / (
                                         3 * (k1 + h ** 2)) - h ** 2 / (3 * (k1 + h ** 2))
                    T[1, 2, 2] = -k1 * h * (2 * np.sqrt(k1) * np.sinh(2 * np.sqrt(k1) * L) + np.sqrt(k1 + h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2))) / (10 * k1 + 2 * h ** 2)
                    T[1, 3, 3] = -h * np.sin(L * np.sqrt(k1 + h ** 2)) / (2 * np.sqrt(k1 + h ** 2))
                    T[1, 4, 4] = -h * (6 * L * h ** 2 * (k1 + h ** 2) ** 8 * (2 * k1 + h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 3 * L * (k1 + h ** 2) ** 9 * (k1 + 2 * h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 6 * h ** 2 * (k1 + h ** 2) ** (15 / 2) * (2 * k1 + h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) - h ** 2 * (k1 + h ** 2) ** (15 / 2) * (
                                               72 * k1 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) - 72 * k1 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) ** 5 + 144 * k1 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) ** 3 - 64 * k1 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) + 16 * k1 + 36 * h ** 2 * np.sin(
                                           L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) - 36 * h ** 2 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) ** 5 + 72 * h ** 2 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) ** 3 - 32 * h ** 2 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) + 8 * h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) + 6 * (
                                               k1 + h ** 2) ** (19 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) + 3 * (
                                               k1 + h ** 2) ** (17 / 2) * (k1 + 2 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) - 2 * (k1 + h ** 2) ** (17 / 2) * (
                                               3 * k1 - h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) + 7 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2))) / (6 * (k1 + h ** 2) ** 10)
                    T[2, 0, 2] = np.sqrt(k1) * h * np.sin(L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) / np.sqrt(k1 + h ** 2)
                    T[2, 0, 3] = h * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) / np.sqrt(k1 + h ** 2) - h * np.sinh(
                        np.sqrt(k1) * L) / np.sqrt(k1)
                    T[2, 1, 2] = -np.sqrt(k1) * h * (np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.sinh(np.sqrt(k1) * L) / (k1 + h ** 2)
                    T[2, 1, 3] = -h * (np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.cosh(np.sqrt(k1) * L) / (k1 + h ** 2)
                    T[2, 2, 4] = -np.sqrt(k1) * (-2 * L * h ** 2 * (k1 + h ** 2) ** (3 / 2) + L * (k1 + h ** 2) ** (
                            5 / 2) + 2 * h ** 2 * (k1 + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2))) * np.sinh(
                        np.sqrt(k1) * L) / (2 * (k1 + h ** 2) ** (5 / 2))
                    T[2, 3, 4] = -(1048576 * k1 ** (25 / 2) * L * (k1 + h ** 2) ** 3 * np.cosh(
                        np.sqrt(k1) * L) - 1048576 * k1 ** (25 / 2) * h ** 2 * (k1 + h ** 2) ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 1048576 * k1 ** (25 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 2 * np.cosh(np.sqrt(k1) * L) + 3932160 * k1 ** (23 / 2) * L * (
                                           k1 + h ** 2) ** 4 * np.cosh(np.sqrt(k1) * L) + 524288 * k1 ** (
                                           23 / 2) * h ** 2 * (k1 + h ** 2) ** (5 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) - 3670016 * k1 ** (23 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 3 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 3670016 * k1 ** (23 / 2) * h ** 2 * (k1 + h ** 2) ** 3 * np.cosh(
                        np.sqrt(k1) * L) + 6684672 * k1 ** (21 / 2) * L * (k1 + h ** 2) ** 5 * np.cosh(
                        np.sqrt(k1) * L) + 1835008 * k1 ** (21 / 2) * h ** 2 * (k1 + h ** 2) ** (7 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) - 5767168 * k1 ** (21 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 4 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 5767168 * k1 ** (21 / 2) * h ** 2 * (k1 + h ** 2) ** 4 * np.cosh(
                        np.sqrt(k1) * L) + 6815744 * k1 ** (19 / 2) * L * (k1 + h ** 2) ** 6 * np.cosh(
                        np.sqrt(k1) * L) + 2883584 * k1 ** (19 / 2) * h ** 2 * (k1 + h ** 2) ** (9 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) - 5373952 * k1 ** (19 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 5 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 5373952 * k1 ** (19 / 2) * h ** 2 * (k1 + h ** 2) ** 5 * np.cosh(
                        np.sqrt(k1) * L) + 4644864 * k1 ** (17 / 2) * L * (k1 + h ** 2) ** 7 * np.cosh(
                        np.sqrt(k1) * L) + 2686976 * k1 ** (17 / 2) * h ** 2 * (k1 + h ** 2) ** (11 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) - 3301376 * k1 ** (17 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 6 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 3301376 * k1 ** (17 / 2) * h ** 2 * (k1 + h ** 2) ** 6 * np.cosh(
                        np.sqrt(k1) * L) + 2230272 * k1 ** (15 / 2) * L * (k1 + h ** 2) ** 8 * np.cosh(
                        np.sqrt(k1) * L) + 1650688 * k1 ** (15 / 2) * h ** 2 * (k1 + h ** 2) ** (13 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) - 1404928 * k1 ** (15 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 7 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 1404928 * k1 ** (15 / 2) * h ** 2 * (k1 + h ** 2) ** 7 * np.cosh(
                        np.sqrt(k1) * L) + 774144 * k1 ** (13 / 2) * L * (k1 + h ** 2) ** 9 * np.cosh(
                        np.sqrt(k1) * L) + 702464 * k1 ** (13 / 2) * h ** 2 * (k1 + h ** 2) ** (15 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) - 422912 * k1 ** (13 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 8 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 422912 * k1 ** (13 / 2) * h ** 2 * (k1 + h ** 2) ** 8 * np.cosh(
                        np.sqrt(k1) * L) + 195840 * k1 ** (11 / 2) * L * (k1 + h ** 2) ** 10 * np.cosh(
                        np.sqrt(k1) * L) + 211456 * k1 ** (11 / 2) * h ** 2 * (k1 + h ** 2) ** (17 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) - 90112 * k1 ** (11 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 9 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 90112 * k1 ** (11 / 2) * h ** 2 * (k1 + h ** 2) ** 9 * np.cosh(
                        np.sqrt(k1) * L) + 35856 * k1 ** (9 / 2) * L * (k1 + h ** 2) ** 11 * np.cosh(
                        np.sqrt(k1) * L) + 45056 * k1 ** (9 / 2) * h ** 2 * (k1 + h ** 2) ** (19 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) - 13328 * k1 ** (9 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 10 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 13328 * k1 ** (9 / 2) * h ** 2 * (k1 + h ** 2) ** 10 * np.cosh(
                        np.sqrt(k1) * L) + 4636 * k1 ** (7 / 2) * L * (k1 + h ** 2) ** 12 * np.cosh(
                        np.sqrt(k1) * L) + 6664 * k1 ** (7 / 2) * h ** 2 * (k1 + h ** 2) ** (21 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) - 1304 * k1 ** (7 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 11 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 1304 * k1 ** (7 / 2) * h ** 2 * (k1 + h ** 2) ** 11 * np.cosh(
                        np.sqrt(k1) * L) + 402 * k1 ** (5 / 2) * L * (k1 + h ** 2) ** 13 * np.cosh(
                        np.sqrt(k1) * L) + 652 * k1 ** (5 / 2) * h ** 2 * (k1 + h ** 2) ** (23 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) - 76 * k1 ** (5 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 12 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 76 * k1 ** (5 / 2) * h ** 2 * (k1 + h ** 2) ** 12 * np.cosh(
                        np.sqrt(k1) * L) + 21 * k1 ** (3 / 2) * L * (k1 + h ** 2) ** 14 * np.cosh(np.sqrt(k1) * L) + 38 * k1 ** (
                                           3 / 2) * h ** 2 * (k1 + h ** 2) ** (25 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) - 2 * k1 ** (3 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 13 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 2 * k1 ** (3 / 2) * h ** 2 * (k1 + h ** 2) ** 13 * np.cosh(np.sqrt(k1) * L) + np.sqrt(
                        k1) * L * (k1 + h ** 2) ** 15 * np.cosh(np.sqrt(k1) * L) / 2 + np.sqrt(k1) * h ** 2 * (k1 + h ** 2) ** (
                                           27 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 2097152 * k1 ** 13 * L * h ** 2 * (k1 + h ** 2) ** 2 * np.sinh(
                        np.sqrt(k1) * L) + 2097152 * k1 ** 13 * h ** 2 * (k1 + h ** 2) ** (3 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) - 7864320 * k1 ** 12 * L * h ** 2 * (
                                           k1 + h ** 2) ** 3 * np.sinh(np.sqrt(k1) * L) + 7340032 * k1 ** 12 * h ** 2 * (
                                           k1 + h ** 2) ** (5 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) + 1048576 * k1 ** 12 * h ** 2 * (k1 + h ** 2) ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) - 1048576 * k1 ** 12 * h ** 2 * (
                                           k1 + h ** 2) ** 2 * np.sinh(np.sqrt(k1) * L) - 1048576 * k1 ** 12 * (
                                           k1 + h ** 2) ** 3 * np.sinh(
                        np.sqrt(k1) * L) - 13369344 * k1 ** 11 * L * h ** 2 * (k1 + h ** 2) ** 4 * np.sinh(
                        np.sqrt(k1) * L) + 11534336 * k1 ** 11 * h ** 2 * (k1 + h ** 2) ** (7 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 3670016 * k1 ** 11 * h ** 2 * (
                                           k1 + h ** 2) ** 3 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 4194304 * k1 ** 11 * h ** 2 * (k1 + h ** 2) ** 3 * np.sinh(
                        np.sqrt(k1) * L) - 3932160 * k1 ** 11 * (k1 + h ** 2) ** 4 * np.sinh(
                        np.sqrt(k1) * L) - 13631488 * k1 ** 10 * L * h ** 2 * (k1 + h ** 2) ** 5 * np.sinh(
                        np.sqrt(k1) * L) + 10747904 * k1 ** 10 * h ** 2 * (k1 + h ** 2) ** (9 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 5767168 * k1 ** 10 * h ** 2 * (
                                           k1 + h ** 2) ** 4 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 7602176 * k1 ** 10 * h ** 2 * (k1 + h ** 2) ** 4 * np.sinh(
                        np.sqrt(k1) * L) - 6684672 * k1 ** 10 * (k1 + h ** 2) ** 5 * np.sinh(
                        np.sqrt(k1) * L) - 9289728 * k1 ** 9 * L * h ** 2 * (k1 + h ** 2) ** 6 * np.sinh(
                        np.sqrt(k1) * L) + 6602752 * k1 ** 9 * h ** 2 * (k1 + h ** 2) ** (11 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 5373952 * k1 ** 9 * h ** 2 * (
                                           k1 + h ** 2) ** 5 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 8257536 * k1 ** 9 * h ** 2 * (k1 + h ** 2) ** 5 * np.sinh(
                        np.sqrt(k1) * L) - 6815744 * k1 ** 9 * (k1 + h ** 2) ** 6 * np.sinh(
                        np.sqrt(k1) * L) - 4460544 * k1 ** 8 * L * h ** 2 * (k1 + h ** 2) ** 7 * np.sinh(
                        np.sqrt(k1) * L) + 2809856 * k1 ** 8 * h ** 2 * (k1 + h ** 2) ** (13 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 3301376 * k1 ** 8 * h ** 2 * (
                                           k1 + h ** 2) ** 6 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 5988352 * k1 ** 8 * h ** 2 * (k1 + h ** 2) ** 6 * np.sinh(
                        np.sqrt(k1) * L) - 4644864 * k1 ** 8 * (k1 + h ** 2) ** 7 * np.sinh(
                        np.sqrt(k1) * L) - 1548288 * k1 ** 7 * L * h ** 2 * (k1 + h ** 2) ** 8 * np.sinh(
                        np.sqrt(k1) * L) + 845824 * k1 ** 7 * h ** 2 * (k1 + h ** 2) ** (15 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 1404928 * k1 ** 7 * h ** 2 * (
                                           k1 + h ** 2) ** 7 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 3055616 * k1 ** 7 * h ** 2 * (k1 + h ** 2) ** 7 * np.sinh(
                        np.sqrt(k1) * L) - 2230272 * k1 ** 7 * (k1 + h ** 2) ** 8 * np.sinh(
                        np.sqrt(k1) * L) - 391680 * k1 ** 6 * L * h ** 2 * (k1 + h ** 2) ** 9 * np.sinh(
                        np.sqrt(k1) * L) + 180224 * k1 ** 6 * h ** 2 * (k1 + h ** 2) ** (17 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 422912 * k1 ** 6 * h ** 2 * (
                                           k1 + h ** 2) ** 8 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 1125376 * k1 ** 6 * h ** 2 * (k1 + h ** 2) ** 8 * np.sinh(
                        np.sqrt(k1) * L) - 774144 * k1 ** 6 * (k1 + h ** 2) ** 9 * np.sinh(
                        np.sqrt(k1) * L) - 71712 * k1 ** 5 * L * h ** 2 * (k1 + h ** 2) ** 10 * np.sinh(
                        np.sqrt(k1) * L) + 26656 * k1 ** 5 * h ** 2 * (k1 + h ** 2) ** (19 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 90112 * k1 ** 5 * h ** 2 * (
                                           k1 + h ** 2) ** 9 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 301568 * k1 ** 5 * h ** 2 * (k1 + h ** 2) ** 9 * np.sinh(
                        np.sqrt(k1) * L) - 195840 * k1 ** 5 * (k1 + h ** 2) ** 10 * np.sinh(
                        np.sqrt(k1) * L) - 9272 * k1 ** 4 * L * h ** 2 * (k1 + h ** 2) ** 11 * np.sinh(
                        np.sqrt(k1) * L) + 2608 * k1 ** 4 * h ** 2 * (k1 + h ** 2) ** (21 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 13328 * k1 ** 4 * h ** 2 * (
                                           k1 + h ** 2) ** 10 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 58384 * k1 ** 4 * h ** 2 * (k1 + h ** 2) ** 10 * np.sinh(
                        np.sqrt(k1) * L) - 35856 * k1 ** 4 * (k1 + h ** 2) ** 11 * np.sinh(
                        np.sqrt(k1) * L) - 804 * k1 ** 3 * L * h ** 2 * (k1 + h ** 2) ** 12 * np.sinh(
                        np.sqrt(k1) * L) + 152 * k1 ** 3 * h ** 2 * (k1 + h ** 2) ** (23 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 1304 * k1 ** 3 * h ** 2 * (
                                           k1 + h ** 2) ** 11 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 7968 * k1 ** 3 * h ** 2 * (k1 + h ** 2) ** 11 * np.sinh(
                        np.sqrt(k1) * L) - 4636 * k1 ** 3 * (k1 + h ** 2) ** 12 * np.sinh(
                        np.sqrt(k1) * L) - 42 * k1 ** 2 * L * h ** 2 * (k1 + h ** 2) ** 13 * np.sinh(
                        np.sqrt(k1) * L) + 4 * k1 ** 2 * h ** 2 * (k1 + h ** 2) ** (25 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 76 * k1 ** 2 * h ** 2 * (k1 + h ** 2) ** 12 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) - 728 * k1 ** 2 * h ** 2 * (
                                           k1 + h ** 2) ** 12 * np.sinh(np.sqrt(k1) * L) - 402 * k1 ** 2 * (
                                           k1 + h ** 2) ** 13 * np.sinh(np.sqrt(k1) * L) - k1 * L * h ** 2 * (
                                           k1 + h ** 2) ** 14 * np.sinh(np.sqrt(k1) * L) + 2 * k1 * h ** 2 * (
                                           k1 + h ** 2) ** 13 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 40 * k1 * h ** 2 * (k1 + h ** 2) ** 13 * np.sinh(np.sqrt(k1) * L) - 21 * k1 * (
                                           k1 + h ** 2) ** 14 * np.sinh(np.sqrt(k1) * L) - h ** 2 * (
                                           k1 + h ** 2) ** 14 * np.sinh(np.sqrt(k1) * L) - (k1 + h ** 2) ** 15 * np.sinh(
                        np.sqrt(k1) * L) / 2) / (np.sqrt(k1) * (k1 + h ** 2) ** 3 * (
                            2097152 * k1 ** 12 + 7864320 * k1 ** 11 * (k1 + h ** 2) + 13369344 * k1 ** 10 * (
                            k1 + h ** 2) ** 2 + 13631488 * k1 ** 9 * (k1 + h ** 2) ** 3 + 9289728 * k1 ** 8 * (
                                    k1 + h ** 2) ** 4 + 4460544 * k1 ** 7 * (
                                    k1 + h ** 2) ** 5 + 1548288 * k1 ** 6 * (
                                    k1 + h ** 2) ** 6 + 391680 * k1 ** 5 * (
                                    k1 + h ** 2) ** 7 + 71712 * k1 ** 4 * (
                                    k1 + h ** 2) ** 8 + 9272 * k1 ** 3 * (k1 + h ** 2) ** 9 + 804 * k1 ** 2 * (
                                    k1 + h ** 2) ** 10 + 42 * k1 * (k1 + h ** 2) ** 11 + (k1 + h ** 2) ** 12))
                    T[3, 0, 2] = np.sqrt(k1) * h * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + k1 * h * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) / np.sqrt(k1 + h ** 2)
                    T[3, 0, 3] = h * (np.sqrt(k1) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + np.sqrt(k1 + h ** 2) * (
                            np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.cosh(np.sqrt(k1) * L)) / np.sqrt(k1 + h ** 2)
                    T[3, 1, 2] = h * (
                            np.sqrt(k1) * (k1 + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) - k1 * np.sqrt(
                        k1 + h ** 2) * (np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.cosh(np.sqrt(k1) * L)) / (k1 + h ** 2) ** (
                                         3 / 2)
                    T[3, 1, 3] = -h * (
                            np.sqrt(k1) * np.sqrt(k1 + h ** 2) * (np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.sinh(np.sqrt(k1) * L) - (
                            k1 + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L)) / (k1 + h ** 2) ** (
                                         3 / 2)
                    T[3, 2, 4] = -(np.sqrt(k1) * (k1 + h ** 2) ** (3 / 2) * (
                            k1 + 2 * h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) - h ** 2) * np.sinh(np.sqrt(k1) * L) + k1 * (
                                           -2 * L * h ** 2 * (k1 + h ** 2) ** (3 / 2) + L * (k1 + h ** 2) ** (
                                           5 / 2) + 2 * h ** 2 * (k1 + h ** 2) * np.sin(
                                       L * np.sqrt(k1 + h ** 2))) * np.cosh(np.sqrt(k1) * L)) / (
                                         2 * (k1 + h ** 2) ** (5 / 2))
                    T[3, 3, 4] = -(1048576 * k1 ** (25 / 2) * L * (k1 + h ** 2) ** 3 * np.sinh(
                        np.sqrt(k1) * L) + 1048576 * k1 ** (25 / 2) * h ** 2 * (k1 + h ** 2) ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) - 1048576 * k1 ** (25 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 2 * np.sinh(np.sqrt(k1) * L) + 3932160 * k1 ** (23 / 2) * L * (
                                           k1 + h ** 2) ** 4 * np.sinh(np.sqrt(k1) * L) - 524288 * k1 ** (
                                           23 / 2) * h ** 2 * (k1 + h ** 2) ** (5 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 3670016 * k1 ** (23 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 3 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 4194304 * k1 ** (23 / 2) * h ** 2 * (k1 + h ** 2) ** 3 * np.sinh(
                        np.sqrt(k1) * L) + 6684672 * k1 ** (21 / 2) * L * (k1 + h ** 2) ** 5 * np.sinh(
                        np.sqrt(k1) * L) - 1835008 * k1 ** (21 / 2) * h ** 2 * (k1 + h ** 2) ** (7 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 5767168 * k1 ** (21 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 4 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 7602176 * k1 ** (21 / 2) * h ** 2 * (k1 + h ** 2) ** 4 * np.sinh(
                        np.sqrt(k1) * L) + 6815744 * k1 ** (19 / 2) * L * (k1 + h ** 2) ** 6 * np.sinh(
                        np.sqrt(k1) * L) - 2883584 * k1 ** (19 / 2) * h ** 2 * (k1 + h ** 2) ** (9 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 5373952 * k1 ** (19 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 5 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 8257536 * k1 ** (19 / 2) * h ** 2 * (k1 + h ** 2) ** 5 * np.sinh(
                        np.sqrt(k1) * L) + 4644864 * k1 ** (17 / 2) * L * (k1 + h ** 2) ** 7 * np.sinh(
                        np.sqrt(k1) * L) - 2686976 * k1 ** (17 / 2) * h ** 2 * (k1 + h ** 2) ** (11 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 3301376 * k1 ** (17 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 6 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 5988352 * k1 ** (17 / 2) * h ** 2 * (k1 + h ** 2) ** 6 * np.sinh(
                        np.sqrt(k1) * L) + 2230272 * k1 ** (15 / 2) * L * (k1 + h ** 2) ** 8 * np.sinh(
                        np.sqrt(k1) * L) - 1650688 * k1 ** (15 / 2) * h ** 2 * (k1 + h ** 2) ** (13 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 1404928 * k1 ** (15 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 7 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 3055616 * k1 ** (15 / 2) * h ** 2 * (k1 + h ** 2) ** 7 * np.sinh(
                        np.sqrt(k1) * L) + 774144 * k1 ** (13 / 2) * L * (k1 + h ** 2) ** 9 * np.sinh(
                        np.sqrt(k1) * L) - 702464 * k1 ** (13 / 2) * h ** 2 * (k1 + h ** 2) ** (15 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 422912 * k1 ** (13 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 8 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 1125376 * k1 ** (13 / 2) * h ** 2 * (k1 + h ** 2) ** 8 * np.sinh(
                        np.sqrt(k1) * L) + 195840 * k1 ** (11 / 2) * L * (k1 + h ** 2) ** 10 * np.sinh(
                        np.sqrt(k1) * L) - 211456 * k1 ** (11 / 2) * h ** 2 * (k1 + h ** 2) ** (17 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 90112 * k1 ** (11 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 9 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 301568 * k1 ** (11 / 2) * h ** 2 * (k1 + h ** 2) ** 9 * np.sinh(
                        np.sqrt(k1) * L) + 35856 * k1 ** (9 / 2) * L * (k1 + h ** 2) ** 11 * np.sinh(
                        np.sqrt(k1) * L) - 45056 * k1 ** (9 / 2) * h ** 2 * (k1 + h ** 2) ** (19 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 13328 * k1 ** (9 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 10 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 58384 * k1 ** (9 / 2) * h ** 2 * (k1 + h ** 2) ** 10 * np.sinh(
                        np.sqrt(k1) * L) + 4636 * k1 ** (7 / 2) * L * (k1 + h ** 2) ** 12 * np.sinh(
                        np.sqrt(k1) * L) - 6664 * k1 ** (7 / 2) * h ** 2 * (k1 + h ** 2) ** (21 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 1304 * k1 ** (7 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 11 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 7968 * k1 ** (7 / 2) * h ** 2 * (k1 + h ** 2) ** 11 * np.sinh(
                        np.sqrt(k1) * L) + 402 * k1 ** (5 / 2) * L * (k1 + h ** 2) ** 13 * np.sinh(
                        np.sqrt(k1) * L) - 652 * k1 ** (5 / 2) * h ** 2 * (k1 + h ** 2) ** (23 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 76 * k1 ** (5 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 12 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 728 * k1 ** (5 / 2) * h ** 2 * (k1 + h ** 2) ** 12 * np.sinh(
                        np.sqrt(k1) * L) + 21 * k1 ** (3 / 2) * L * (k1 + h ** 2) ** 14 * np.sinh(np.sqrt(k1) * L) - 38 * k1 ** (
                                           3 / 2) * h ** 2 * (k1 + h ** 2) ** (25 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) + 2 * k1 ** (3 / 2) * h ** 2 * (
                                           k1 + h ** 2) ** 13 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.sinh(
                        np.sqrt(k1) * L) - 40 * k1 ** (3 / 2) * h ** 2 * (k1 + h ** 2) ** 13 * np.sinh(np.sqrt(k1) * L) + np.sqrt(
                        k1) * L * (k1 + h ** 2) ** 15 * np.sinh(np.sqrt(k1) * L) / 2 - np.sqrt(k1) * h ** 2 * (k1 + h ** 2) ** (
                                           27 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.sinh(np.sqrt(k1) * L) - np.sqrt(
                        k1) * h ** 2 * (k1 + h ** 2) ** 14 * np.sinh(np.sqrt(k1) * L) - 2097152 * k1 ** 13 * L * h ** 2 * (
                                           k1 + h ** 2) ** 2 * np.cosh(np.sqrt(k1) * L) + 2097152 * k1 ** 13 * h ** 2 * (
                                           k1 + h ** 2) ** (3 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 7864320 * k1 ** 12 * L * h ** 2 * (k1 + h ** 2) ** 3 * np.cosh(
                        np.sqrt(k1) * L) + 8388608 * k1 ** 12 * h ** 2 * (k1 + h ** 2) ** (5 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 1048576 * k1 ** 12 * h ** 2 * (
                                           k1 + h ** 2) ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 1048576 * k1 ** 12 * h ** 2 * (k1 + h ** 2) ** 2 * np.cosh(
                        np.sqrt(k1) * L) - 13369344 * k1 ** 11 * L * h ** 2 * (k1 + h ** 2) ** 4 * np.cosh(
                        np.sqrt(k1) * L) + 15204352 * k1 ** 11 * h ** 2 * (k1 + h ** 2) ** (7 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 4194304 * k1 ** 11 * h ** 2 * (
                                           k1 + h ** 2) ** 3 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 4194304 * k1 ** 11 * h ** 2 * (k1 + h ** 2) ** 3 * np.cosh(
                        np.sqrt(k1) * L) - 13631488 * k1 ** 10 * L * h ** 2 * (k1 + h ** 2) ** 5 * np.cosh(
                        np.sqrt(k1) * L) + 16515072 * k1 ** 10 * h ** 2 * (k1 + h ** 2) ** (9 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 7602176 * k1 ** 10 * h ** 2 * (
                                           k1 + h ** 2) ** 4 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 7602176 * k1 ** 10 * h ** 2 * (k1 + h ** 2) ** 4 * np.cosh(
                        np.sqrt(k1) * L) - 9289728 * k1 ** 9 * L * h ** 2 * (k1 + h ** 2) ** 6 * np.cosh(
                        np.sqrt(k1) * L) + 11976704 * k1 ** 9 * h ** 2 * (k1 + h ** 2) ** (11 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 8257536 * k1 ** 9 * h ** 2 * (
                                           k1 + h ** 2) ** 5 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 8257536 * k1 ** 9 * h ** 2 * (k1 + h ** 2) ** 5 * np.cosh(
                        np.sqrt(k1) * L) - 4460544 * k1 ** 8 * L * h ** 2 * (k1 + h ** 2) ** 7 * np.cosh(
                        np.sqrt(k1) * L) + 6111232 * k1 ** 8 * h ** 2 * (k1 + h ** 2) ** (13 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 5988352 * k1 ** 8 * h ** 2 * (
                                           k1 + h ** 2) ** 6 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 5988352 * k1 ** 8 * h ** 2 * (k1 + h ** 2) ** 6 * np.cosh(
                        np.sqrt(k1) * L) - 1548288 * k1 ** 7 * L * h ** 2 * (k1 + h ** 2) ** 8 * np.cosh(
                        np.sqrt(k1) * L) + 2250752 * k1 ** 7 * h ** 2 * (k1 + h ** 2) ** (15 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 3055616 * k1 ** 7 * h ** 2 * (
                                           k1 + h ** 2) ** 7 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 3055616 * k1 ** 7 * h ** 2 * (k1 + h ** 2) ** 7 * np.cosh(
                        np.sqrt(k1) * L) - 391680 * k1 ** 6 * L * h ** 2 * (k1 + h ** 2) ** 9 * np.cosh(
                        np.sqrt(k1) * L) + 603136 * k1 ** 6 * h ** 2 * (k1 + h ** 2) ** (17 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 1125376 * k1 ** 6 * h ** 2 * (
                                           k1 + h ** 2) ** 8 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 1125376 * k1 ** 6 * h ** 2 * (k1 + h ** 2) ** 8 * np.cosh(
                        np.sqrt(k1) * L) - 71712 * k1 ** 5 * L * h ** 2 * (k1 + h ** 2) ** 10 * np.cosh(
                        np.sqrt(k1) * L) + 116768 * k1 ** 5 * h ** 2 * (k1 + h ** 2) ** (19 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 301568 * k1 ** 5 * h ** 2 * (
                                           k1 + h ** 2) ** 9 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 301568 * k1 ** 5 * h ** 2 * (k1 + h ** 2) ** 9 * np.cosh(
                        np.sqrt(k1) * L) - 9272 * k1 ** 4 * L * h ** 2 * (k1 + h ** 2) ** 11 * np.cosh(
                        np.sqrt(k1) * L) + 15936 * k1 ** 4 * h ** 2 * (k1 + h ** 2) ** (21 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 58384 * k1 ** 4 * h ** 2 * (
                                           k1 + h ** 2) ** 10 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 58384 * k1 ** 4 * h ** 2 * (k1 + h ** 2) ** 10 * np.cosh(
                        np.sqrt(k1) * L) - 804 * k1 ** 3 * L * h ** 2 * (k1 + h ** 2) ** 12 * np.cosh(
                        np.sqrt(k1) * L) + 1456 * k1 ** 3 * h ** 2 * (k1 + h ** 2) ** (23 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 7968 * k1 ** 3 * h ** 2 * (
                                           k1 + h ** 2) ** 11 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 7968 * k1 ** 3 * h ** 2 * (k1 + h ** 2) ** 11 * np.cosh(
                        np.sqrt(k1) * L) - 42 * k1 ** 2 * L * h ** 2 * (k1 + h ** 2) ** 13 * np.cosh(
                        np.sqrt(k1) * L) + 80 * k1 ** 2 * h ** 2 * (k1 + h ** 2) ** (25 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cosh(np.sqrt(k1) * L) + 728 * k1 ** 2 * h ** 2 * (
                                           k1 + h ** 2) ** 12 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 728 * k1 ** 2 * h ** 2 * (k1 + h ** 2) ** 12 * np.cosh(
                        np.sqrt(k1) * L) - k1 * L * h ** 2 * (k1 + h ** 2) ** 14 * np.cosh(np.sqrt(k1) * L) + 2 * k1 * h ** 2 * (
                                           k1 + h ** 2) ** (27 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) + 40 * k1 * h ** 2 * (k1 + h ** 2) ** 13 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - 40 * k1 * h ** 2 * (k1 + h ** 2) ** 13 * np.cosh(np.sqrt(k1) * L) + h ** 2 * (
                                           k1 + h ** 2) ** 14 * np.cos(L * np.sqrt(k1 + h ** 2)) * np.cosh(
                        np.sqrt(k1) * L) - h ** 2 * (k1 + h ** 2) ** 14 * np.cosh(np.sqrt(k1) * L)) / ((k1 + h ** 2) ** 3 * (
                            2097152 * k1 ** 12 + 7864320 * k1 ** 11 * (k1 + h ** 2) + 13369344 * k1 ** 10 * (
                            k1 + h ** 2) ** 2 + 13631488 * k1 ** 9 * (k1 + h ** 2) ** 3 + 9289728 * k1 ** 8 * (
                                    k1 + h ** 2) ** 4 + 4460544 * k1 ** 7 * (
                                    k1 + h ** 2) ** 5 + 1548288 * k1 ** 6 * (
                                    k1 + h ** 2) ** 6 + 391680 * k1 ** 5 * (
                                    k1 + h ** 2) ** 7 + 71712 * k1 ** 4 * (
                                    k1 + h ** 2) ** 8 + 9272 * k1 ** 3 * (k1 + h ** 2) ** 9 + 804 * k1 ** 2 * (
                                    k1 + h ** 2) ** 10 + 42 * k1 * (k1 + h ** 2) ** 11 + (k1 + h ** 2) ** 12))
                if k1 < 0:
                    T[0, 0, 0] = h * (4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) ** 2 + 4 * k1 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 7 * k1 + 2 * h ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) ** 2 + 2 * h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) - 3 * h ** 2 + (
                                              k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2)) ** 2 - 2 * (
                                              k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2))) / (6 * (k1 + h ** 2))
                    T[0, 0, 1] = h * (4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) - 5 * k1 + 2 * h ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 3 * h ** 2 + (k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2))) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) / (3 * (k1 + h ** 2) ** (3 / 2))
                    T[0, 0, 4] = (-6 * L * h ** 2 * (k1 + h ** 2) ** (7 / 2) * (2 * k1 + h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + 3 * L * (k1 + h ** 2) ** (9 / 2) * (k1 + 2 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + 2 * h ** 2 * (k1 + h ** 2) ** 4 * (
                                          np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 + 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) - 2) + h ** 2 * (k1 + h ** 2) ** 3 * (
                                          -6 * k1 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) + 8 * k1 * np.sin(
                                      L * np.sqrt(k1 + h ** 2)) ** 2 + 6 * k1 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) ** 5 - 12 * k1 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) ** 3 - 2 * k1 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) + 8 * k1 - 3 * h ** 2 * np.sin(
                                      L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(L * np.sqrt(k1 + h ** 2)) + 4 * h ** 2 * np.sin(
                                      L * np.sqrt(k1 + h ** 2)) ** 2 + 3 * h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) ** 5 - 6 * h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) ** 3 - h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) + 4 * h ** 2)) / (6 * (k1 + h ** 2) ** 5)
                    T[0, 1, 1] = h * (4 * k1 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 + 8 * k1 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 7 * k1 + 2 * h ** 2 * np.sin(
                        L * np.sqrt(k1 + h ** 2)) ** 2 + 4 * h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) - 3 * h ** 2 + (
                                              k1 + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 - (k1 + h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2))) / (6 * (k1 + h ** 2) ** 2)
                    T[0, 1, 4] = (-L * h ** 2 * (k1 + h ** 2) ** 5 * (2 * k1 + h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - L * (k1 + h ** 2) ** 6 * (k1 + 2 * h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) / 2 + h ** 2 * (k1 + h ** 2) ** (9 / 2) * (
                                          4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) + 2 * k1 + 2 * h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) / 3 + (
                                          k1 + h ** 2) ** (11 / 2) * (
                                          3 * k1 - 2 * h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) + 8 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) / 6) / (k1 + h ** 2) ** 7

                    if 2 * k1 == -k1 - h ** 2:
                        T[0, 2, 2] = np.real(h * (np.cos(L * np.sqrt(k1 + h ** 2)) - np.cos(np.sqrt(2) * L * np.sqrt(k1 + h ** 2))) / 4)

                    elif 4 * k1 == -k1 - h ** 2:
                        T[0, 2, 2] = np.real(L * h * np.sqrt(k1 + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) / 16)

                    else:
                        T[0, 2, 2] = np.real(k1 * h * (np.cos(L * np.sqrt(k1 + h ** 2)) - np.cosh(2 * np.sqrt(k1) * L)) /
                                             (2 * (5 * k1 + h ** 2)))
                    if 4 * k1 != -k1 - h ** 2:
                        T[0, 2, 3] = np.real(-h * (1j * k1 ** (3 / 2) - (-k1) ** (3 / 2)) * np.sin(L * np.sqrt(k1 + h ** 2)) /
                                             np.sqrt(-k1) * np.sqrt(k1 + h ** 2) * (5 * k1 + h ** 2))

                    else:
                        T[0, 2, 3] = 0
                    T[0, 3, 3] = h * (np.cos(L * np.sqrt(k1 + h ** 2)) - 1) / (2 * (k1 + h ** 2))
                    T[0, 4, 4] = h * (-6 * L * h ** 2 * (k1 + h ** 2) ** (15 / 2) * (2 * k1 + h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) - 3 * L * (k1 + h ** 2) ** (17 / 2) * (k1 + 2 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + h ** 2 * (k1 + h ** 2) ** 7 * (
                                              12 * k1 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 6 - 6 * k1 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(L * np.sqrt(k1 + h ** 2)) - 36 * k1 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 4 + 40 * k1 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 2 + 12 * k1 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 6 + 6 * k1 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 5 - 12 * k1 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 3 - 10 * k1 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) + 4 * k1 + 6 * h ** 2 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 6 - 3 * h ** 2 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) - 18 * h ** 2 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 4 + 20 * h ** 2 * np.sin(
                                          L * np.sqrt(k1 + h ** 2)) ** 2 + 6 * h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 6 + 3 * h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 5 - 6 * h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 3 - 5 * h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) + 2 * h ** 2) + 6 * (k1 + h ** 2) ** 9 * (
                                              np.cos(L * np.sqrt(k1 + h ** 2)) - 1) + (k1 + h ** 2) ** 8 * (
                                              -6 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) + 6 * k1 + h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) ** 2 - 14 * h ** 2 * np.cos(
                                          L * np.sqrt(k1 + h ** 2)) + 13 * h ** 2)) / (6 * (k1 + h ** 2) ** 10)
                    T[1, 0, 0] = -h * (
                            4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) + k1 + 2 * h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) + (
                            k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2))) * np.sin(L * np.sqrt(k1 + h ** 2)) / (
                                         3 * np.sqrt(k1 + h ** 2))
                    T[1, 0, 1] = h * (-(5 * k1 + 3 * h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 + (
                            4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) - 5 * k1 + 2 * h ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 3 * h ** 2 + (k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2))) * np.cos(
                        L * np.sqrt(k1 + h ** 2))) / (3 * (k1 + h ** 2))
                    T[1, 0, 4] = (-L * h ** 2 * (k1 + h ** 2) ** 4 * (2 * k1 + h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + L * (k1 + h ** 2) ** 5 * (k1 + 2 * h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) / 2 + 2 * h ** 2 * (k1 + h ** 2) ** (9 / 2) * (
                                          np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) / 3 - h ** 2 * (k1 + h ** 2) ** (7 / 2) * (2 * k1 + h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + 2 * h ** 2 * (k1 + h ** 2) ** (7 / 2) * (
                                          4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) + 2 * k1 + 2 * h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) / 3 + (
                                          k1 + h ** 2) ** (9 / 2) * (k1 + 2 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) / 2) / (k1 + h ** 2) ** 5
                    T[1, 1, 1] = h * (8 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) - 7 * k1 + 4 * h ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 3 * h ** 2 + 2 * (k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2))) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) / (6 * (k1 + h ** 2) ** (3 / 2))
                    T[1, 1, 4] = 2 * k1 * L * h ** 2 * np.sin(L * np.sqrt(k1 + h ** 2)) / (k1 + h ** 2) ** (
                            3 / 2) + k1 * L * np.sin(L * np.sqrt(k1 + h ** 2)) / (
                                         2 * np.sqrt(k1 + h ** 2)) - 8 * k1 * h ** 2 * np.sin(
                        L * np.sqrt(k1 + h ** 2)) ** 2 / (3 * (k1 + h ** 2) ** 2) - 4 * k1 * h ** 2 * np.cos(
                        L * np.sqrt(k1 + h ** 2)) / (3 * (k1 + h ** 2) ** 2) + 4 * k1 * h ** 2 / (
                                         3 * (k1 + h ** 2) ** 2) + L * h ** 4 * np.sin(L * np.sqrt(k1 + h ** 2)) / (
                                         k1 + h ** 2) ** (3 / 2) + L * h ** 2 * np.sin(L * np.sqrt(k1 + h ** 2)) / np.sqrt(
                        k1 + h ** 2) - 4 * h ** 4 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 / (
                                         3 * (k1 + h ** 2) ** 2) - 2 * h ** 4 * np.cos(L * np.sqrt(k1 + h ** 2)) / (
                                         3 * (k1 + h ** 2) ** 2) + 2 * h ** 4 / (
                                         3 * (k1 + h ** 2) ** 2) + 2 * h ** 2 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 2 / (
                                         3 * (k1 + h ** 2)) + h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) / (
                                         3 * (k1 + h ** 2)) - h ** 2 / (3 * (k1 + h ** 2))

                    if 2 * k1 == -k1 - h ** 2:
                        T[1, 2, 2] = np.real(h * np.sqrt(k1 + h ** 2) * (-np.sin(L * np.sqrt(k1 + h ** 2)) + np.sqrt(2) *
                                                                      np.sin(np.sqrt(2) * L * np.sqrt(k1 + h ** 2))))
                    elif 4 * k1 == -k1 - h ** 2:
                        T[1, 2, 2] = np.real(h * (L * (k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2)) + np.sqrt(k1 + h ** 2) *
                                                  np.sin(L * np.sqrt(k1 + h ** 2))))
                    else:
                        T[1, 2, 2] = np.real(-k1 * h * (2 * np.sqrt(k1) * np.sinh(2 * np.sqrt(k1) * L) + np.sqrt(k1 + h ** 2) *
                                                        np.sin(L * np.sqrt(k1 + h ** 2))) / (10 * k1 + 2 * h ** 2))

                    if 4 * k1 != -k1 - h ** 2:
                        T[1, 2, 3] = np.real((-h * (1j * k1 ** (3 / 2) - (-k1) ** (3 / 2)) * np.cos(L * np.sqrt(k1 + h ** 2)) / (
                                    np.sqrt(-k1) * (5 * k1 + h ** 2))))
                    else:
                        T[1, 2, 3] = 0
                    T[1, 3, 3] = -h * np.sin(L * np.sqrt(k1 + h ** 2)) / (2 * np.sqrt(k1 + h ** 2))
                    T[1, 4, 4] = -h * (6 * L * h ** 2 * (k1 + h ** 2) ** 8 * (2 * k1 + h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 3 * L * (k1 + h ** 2) ** 9 * (k1 + 2 * h ** 2) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 6 * h ** 2 * (k1 + h ** 2) ** (15 / 2) * (2 * k1 + h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) - h ** 2 * (k1 + h ** 2) ** (15 / 2) * (
                                               72 * k1 * np.sin(L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) - 72 * k1 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) ** 5 + 144 * k1 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) ** 3 - 64 * k1 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) + 16 * k1 + 36 * h ** 2 * np.sin(
                                           L * np.sqrt(k1 + h ** 2)) ** 4 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) - 36 * h ** 2 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) ** 5 + 72 * h ** 2 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) ** 3 - 32 * h ** 2 * np.cos(
                                           L * np.sqrt(k1 + h ** 2)) + 8 * h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) + 6 * (
                                               k1 + h ** 2) ** (19 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) + 3 * (
                                               k1 + h ** 2) ** (17 / 2) * (k1 + 2 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) - 2 * (k1 + h ** 2) ** (17 / 2) * (
                                               3 * k1 - h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) + 7 * h ** 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2))) / (6 * (k1 + h ** 2) ** 10)
                    T[2, 0, 2] = -h * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) / np.sqrt(k1 + h ** 2)
                    T[2, 0, 3] = h * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) / np.sqrt(k1 + h ** 2) - h * np.sin(
                        L * np.sqrt(-k1)) / np.sqrt(-k1)
                    T[2, 1, 2] = h * np.sqrt(-k1) * (np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.sin(L * np.sqrt(-k1)) / (k1 + h ** 2)
                    T[2, 1, 3] = -h * (np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.cos(L * np.sqrt(-k1)) / (k1 + h ** 2)
                    T[2, 2, 4] = np.sqrt(-k1) * (-2 * L * h ** 2 * (k1 + h ** 2) ** (3 / 2) + L * (k1 + h ** 2) ** (
                            5 / 2) + 2 * h ** 2 * (k1 + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2))) * np.sin(
                        L * np.sqrt(-k1)) / (2 * (k1 + h ** 2) ** (5 / 2))
                    T[2, 3, 4] = -(-33554432 * k1 ** 15 * L * h ** 2 * (k1 + h ** 2) ** 2 * np.sin(
                        L * np.sqrt(-k1)) + 33554432 * k1 ** 15 * h ** 2 * (k1 + h ** 2) ** (3 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 142606336 * k1 ** 14 * L * h ** 2 * (
                                           k1 + h ** 2) ** 3 * np.sin(
                        L * np.sqrt(-k1)) + 134217728 * k1 ** 14 * h ** 2 * (k1 + h ** 2) ** (5 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) + 16777216 * k1 ** 14 * h ** 2 * (
                                           k1 + h ** 2) ** 2 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 16777216 * k1 ** 14 * h ** 2 * (k1 + h ** 2) ** 2 * np.sin(
                        L * np.sqrt(-k1)) - 16777216 * k1 ** 14 * (k1 + h ** 2) ** 3 * np.sin(
                        L * np.sqrt(-k1)) - 278921216 * k1 ** 13 * L * h ** 2 * (k1 + h ** 2) ** 4 * np.sin(
                        L * np.sqrt(-k1)) + 245366784 * k1 ** 13 * h ** 2 * (k1 + h ** 2) ** (7 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) + 67108864 * k1 ** 13 * h ** 2 * (
                                           k1 + h ** 2) ** 3 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 75497472 * k1 ** 13 * h ** 2 * (k1 + h ** 2) ** 3 * np.sin(
                        L * np.sqrt(-k1)) - 71303168 * k1 ** 13 * (k1 + h ** 2) ** 4 * np.sin(
                        L * np.sqrt(-k1)) - 332922880 * k1 ** 12 * L * h ** 2 * (k1 + h ** 2) ** 5 * np.sin(
                        L * np.sqrt(-k1)) + 271581184 * k1 ** 12 * h ** 2 * (k1 + h ** 2) ** (9 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) + 122683392 * k1 ** 12 * h ** 2 * (
                                           k1 + h ** 2) ** 4 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 156237824 * k1 ** 12 * h ** 2 * (k1 + h ** 2) ** 4 * np.sin(
                        L * np.sqrt(-k1)) - 139460608 * k1 ** 12 * (k1 + h ** 2) ** 5 * np.sin(
                        L * np.sqrt(-k1)) - 271056896 * k1 ** 11 * L * h ** 2 * (k1 + h ** 2) ** 6 * np.sin(
                        L * np.sqrt(-k1)) + 203161600 * k1 ** 11 * h ** 2 * (k1 + h ** 2) ** (11 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) + 135790592 * k1 ** 11 * h ** 2 * (
                                           k1 + h ** 2) ** 5 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 197132288 * k1 ** 11 * h ** 2 * (k1 + h ** 2) ** 5 * np.sin(
                        L * np.sqrt(-k1)) - 166461440 * k1 ** 11 * (k1 + h ** 2) ** 6 * np.sin(
                        L * np.sqrt(-k1)) - 159318016 * k1 ** 10 * L * h ** 2 * (k1 + h ** 2) ** 7 * np.sin(
                        L * np.sqrt(-k1)) + 108527616 * k1 ** 10 * h ** 2 * (k1 + h ** 2) ** (13 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) + 101580800 * k1 ** 10 * h ** 2 * (
                                           k1 + h ** 2) ** 6 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 169476096 * k1 ** 10 * h ** 2 * (k1 + h ** 2) ** 6 * np.sin(
                        L * np.sqrt(-k1)) - 135528448 * k1 ** 10 * (k1 + h ** 2) ** 7 * np.sin(
                        L * np.sqrt(-k1)) - 69746688 * k1 ** 9 * L * h ** 2 * (k1 + h ** 2) ** 8 * np.sin(
                        L * np.sqrt(-k1)) + 42614784 * k1 ** 9 * h ** 2 * (k1 + h ** 2) ** (15 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) + 54263808 * k1 ** 9 * h ** 2 * (
                                           k1 + h ** 2) ** 7 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 105054208 * k1 ** 9 * h ** 2 * (k1 + h ** 2) ** 7 * np.sin(
                        L * np.sqrt(-k1)) - 79659008 * k1 ** 9 * (k1 + h ** 2) ** 8 * np.sin(
                        L * np.sqrt(-k1)) - 23113728 * k1 ** 8 * L * h ** 2 * (k1 + h ** 2) ** 9 * np.sin(
                        L * np.sqrt(-k1)) + 12460032 * k1 ** 8 * h ** 2 * (k1 + h ** 2) ** (17 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) + 21307392 * k1 ** 8 * h ** 2 * (
                                           k1 + h ** 2) ** 8 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 48439296 * k1 ** 8 * h ** 2 * (k1 + h ** 2) ** 8 * np.sin(
                        L * np.sqrt(-k1)) - 34873344 * k1 ** 8 * (k1 + h ** 2) ** 9 * np.sin(
                        L * np.sqrt(-k1)) - 5829120 * k1 ** 7 * L * h ** 2 * (k1 + h ** 2) ** 10 * np.sin(
                        L * np.sqrt(-k1)) + 2714112 * k1 ** 7 * h ** 2 * (k1 + h ** 2) ** (19 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) + 6230016 * k1 ** 7 * h ** 2 * (
                                           k1 + h ** 2) ** 9 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 16883712 * k1 ** 7 * h ** 2 * (k1 + h ** 2) ** 9 * np.sin(
                        L * np.sqrt(-k1)) - 11556864 * k1 ** 7 * (k1 + h ** 2) ** 10 * np.sin(
                        L * np.sqrt(-k1)) - 1113728 * k1 ** 6 * L * h ** 2 * (k1 + h ** 2) ** 11 * np.sin(
                        L * np.sqrt(-k1)) + 435200 * k1 ** 6 * h ** 2 * (k1 + h ** 2) ** (21 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) + 1357056 * k1 ** 6 * h ** 2 * (
                                           k1 + h ** 2) ** 10 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 4472064 * k1 ** 6 * h ** 2 * (k1 + h ** 2) ** 10 * np.sin(
                        L * np.sqrt(-k1)) - 2914560 * k1 ** 6 * (k1 + h ** 2) ** 11 * np.sin(
                        L * np.sqrt(-k1)) - 158752 * k1 ** 5 * L * h ** 2 * (k1 + h ** 2) ** 12 * np.sin(
                        L * np.sqrt(-k1)) + 49952 * k1 ** 5 * h ** 2 * (k1 + h ** 2) ** (23 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) + 217600 * k1 ** 5 * h ** 2 * (
                                           k1 + h ** 2) ** 11 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 896128 * k1 ** 5 * h ** 2 * (k1 + h ** 2) ** 11 * np.sin(
                        L * np.sqrt(-k1)) - 556864 * k1 ** 5 * (k1 + h ** 2) ** 12 * np.sin(
                        L * np.sqrt(-k1)) - 16376 * k1 ** 4 * L * h ** 2 * (k1 + h ** 2) ** 13 * np.sin(
                        L * np.sqrt(-k1)) + 3888 * k1 ** 4 * h ** 2 * (k1 + h ** 2) ** (25 / 2) * np.sin(L * np.sqrt(-k1)) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + 24976 * k1 ** 4 * h ** 2 * (k1 + h ** 2) ** 12 * np.sin(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) - 133776 * k1 ** 4 * h ** 2 * (
                                           k1 + h ** 2) ** 12 * np.sin(L * np.sqrt(-k1)) - 79376 * k1 ** 4 * (
                                           k1 + h ** 2) ** 13 * np.sin(L * np.sqrt(-k1)) - 1156 * k1 ** 3 * L * h ** 2 * (

                                           k1 + h ** 2) ** 14 * np.sin(L * np.sqrt(-k1)) + 184 * k1 ** 3 * h ** 2 * (
                                           k1 + h ** 2) ** (27 / 2) * np.sin(L * np.sqrt(-k1)) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + 1944 * k1 ** 3 * h ** 2 * (k1 + h ** 2) ** 13 * np.sin(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) - 14432 * k1 ** 3 * h ** 2 * (
                                           k1 + h ** 2) ** 13 * np.sin(L * np.sqrt(-k1)) - 8188 * k1 ** 3 * (
                                           k1 + h ** 2) ** 14 * np.sin(L * np.sqrt(-k1)) - 50 * k1 ** 2 * L * h ** 2 * (
                                           k1 + h ** 2) ** 15 * np.sin(L * np.sqrt(-k1)) + 4 * k1 ** 2 * h ** 2 * (
                                           k1 + h ** 2) ** (29 / 2) * np.sin(L * np.sqrt(-k1)) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + 92 * k1 ** 2 * h ** 2 * (k1 + h ** 2) ** 14 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 1064 * k1 ** 2 * h ** 2 * (k1 + h ** 2) ** 14 * np.sin(
                        L * np.sqrt(-k1)) - 578 * k1 ** 2 * (k1 + h ** 2) ** 15 * np.sin(L * np.sqrt(-k1)) - k1 * L * h ** 2 * (
                                           k1 + h ** 2) ** 16 * np.sin(L * np.sqrt(-k1)) + 2 * k1 * h ** 2 * (
                                           k1 + h ** 2) ** 15 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 48 * k1 * h ** 2 * (k1 + h ** 2) ** 15 * np.sin(
                        L * np.sqrt(-k1)) - 25 * k1 * (k1 + h ** 2) ** 16 * np.sin(L * np.sqrt(-k1)) + 16777216 * L * (-k1) ** (
                                           29 / 2) * (k1 + h ** 2) ** 3 * np.cos(L * np.sqrt(-k1)) - 71303168 * L * (
                                       -k1) ** (27 / 2) * (k1 + h ** 2) ** 4 * np.cos(L * np.sqrt(-k1)) + 139460608 * L * (
                                       -k1) ** (25 / 2) * (k1 + h ** 2) ** 5 * np.cos(L * np.sqrt(-k1)) - 166461440 * L * (
                                       -k1) ** (23 / 2) * (k1 + h ** 2) ** 6 * np.cos(L * np.sqrt(-k1)) + 135528448 * L * (
                                       -k1) ** (21 / 2) * (k1 + h ** 2) ** 7 * np.cos(L * np.sqrt(-k1)) - 79659008 * L * (
                                       -k1) ** (19 / 2) * (k1 + h ** 2) ** 8 * np.cos(L * np.sqrt(-k1)) + 34873344 * L * (
                                       -k1) ** (17 / 2) * (k1 + h ** 2) ** 9 * np.cos(L * np.sqrt(-k1)) - 11556864 * L * (
                                       -k1) ** (15 / 2) * (k1 + h ** 2) ** 10 * np.cos(L * np.sqrt(-k1)) + 2914560 * L * (
                                       -k1) ** (13 / 2) * (k1 + h ** 2) ** 11 * np.cos(L * np.sqrt(-k1)) - 556864 * L * (
                                       -k1) ** (11 / 2) * (k1 + h ** 2) ** 12 * np.cos(L * np.sqrt(-k1)) + 79376 * L * (
                                       -k1) ** (9 / 2) * (k1 + h ** 2) ** 13 * np.cos(L * np.sqrt(-k1)) - 8188 * L * (
                                       -k1) ** (7 / 2) * (k1 + h ** 2) ** 14 * np.cos(L * np.sqrt(-k1)) + 578 * L * (-k1) ** (
                                           5 / 2) * (k1 + h ** 2) ** 15 * np.cos(L * np.sqrt(-k1)) - 25 * L * (-k1) ** (
                                           3 / 2) * (k1 + h ** 2) ** 16 * np.cos(L * np.sqrt(-k1)) + L * np.sqrt(-k1) * (
                                           k1 + h ** 2) ** 17 * np.cos(L * np.sqrt(-k1)) / 2 - 16777216 * h ** 2 * (
                                       -k1) ** (29 / 2) * (k1 + h ** 2) ** 2 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 16777216 * h ** 2 * (-k1) ** (29 / 2) * (k1 + h ** 2) ** 2 * np.cos(
                        L * np.sqrt(-k1)) - 8388608 * h ** 2 * (-k1) ** (27 / 2) * (k1 + h ** 2) ** (5 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) + 67108864 * h ** 2 * (-k1) ** (27 / 2) * (
                                           k1 + h ** 2) ** 3 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 67108864 * h ** 2 * (-k1) ** (27 / 2) * (k1 + h ** 2) ** 3 * np.cos(
                        L * np.sqrt(-k1)) + 33554432 * h ** 2 * (-k1) ** (25 / 2) * (k1 + h ** 2) ** (7 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) - 122683392 * h ** 2 * (-k1) ** (25 / 2) * (
                                           k1 + h ** 2) ** 4 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 122683392 * h ** 2 * (-k1) ** (25 / 2) * (k1 + h ** 2) ** 4 * np.cos(
                        L * np.sqrt(-k1)) - 61341696 * h ** 2 * (-k1) ** (23 / 2) * (k1 + h ** 2) ** (9 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) + 135790592 * h ** 2 * (-k1) ** (23 / 2) * (
                                           k1 + h ** 2) ** 5 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 135790592 * h ** 2 * (-k1) ** (23 / 2) * (k1 + h ** 2) ** 5 * np.cos(
                        L * np.sqrt(-k1)) + 67895296 * h ** 2 * (-k1) ** (21 / 2) * (k1 + h ** 2) ** (11 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) - 101580800 * h ** 2 * (-k1) ** (21 / 2) * (
                                           k1 + h ** 2) ** 6 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 101580800 * h ** 2 * (-k1) ** (21 / 2) * (k1 + h ** 2) ** 6 * np.cos(
                        L * np.sqrt(-k1)) - 50790400 * h ** 2 * (-k1) ** (19 / 2) * (k1 + h ** 2) ** (13 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) + 54263808 * h ** 2 * (-k1) ** (19 / 2) * (
                                           k1 + h ** 2) ** 7 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 54263808 * h ** 2 * (-k1) ** (19 / 2) * (k1 + h ** 2) ** 7 * np.cos(
                        L * np.sqrt(-k1)) + 27131904 * h ** 2 * (-k1) ** (17 / 2) * (k1 + h ** 2) ** (15 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) - 21307392 * h ** 2 * (-k1) ** (17 / 2) * (
                                           k1 + h ** 2) ** 8 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 21307392 * h ** 2 * (-k1) ** (17 / 2) * (k1 + h ** 2) ** 8 * np.cos(
                        L * np.sqrt(-k1)) - 10653696 * h ** 2 * (-k1) ** (15 / 2) * (k1 + h ** 2) ** (17 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) + 6230016 * h ** 2 * (-k1) ** (15 / 2) * (
                                           k1 + h ** 2) ** 9 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 6230016 * h ** 2 * (-k1) ** (15 / 2) * (k1 + h ** 2) ** 9 * np.cos(
                        L * np.sqrt(-k1)) + 3115008 * h ** 2 * (-k1) ** (13 / 2) * (k1 + h ** 2) ** (19 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) - 1357056 * h ** 2 * (-k1) ** (13 / 2) * (
                                           k1 + h ** 2) ** 10 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 1357056 * h ** 2 * (-k1) ** (13 / 2) * (k1 + h ** 2) ** 10 * np.cos(
                        L * np.sqrt(-k1)) - 678528 * h ** 2 * (-k1) ** (11 / 2) * (k1 + h ** 2) ** (21 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) + 217600 * h ** 2 * (-k1) ** (11 / 2) * (
                                           k1 + h ** 2) ** 11 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 217600 * h ** 2 * (-k1) ** (11 / 2) * (k1 + h ** 2) ** 11 * np.cos(
                        L * np.sqrt(-k1)) + 108800 * h ** 2 * (-k1) ** (9 / 2) * (k1 + h ** 2) ** (23 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) - 24976 * h ** 2 * (-k1) ** (9 / 2) * (
                                           k1 + h ** 2) ** 12 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 24976 * h ** 2 * (-k1) ** (9 / 2) * (k1 + h ** 2) ** 12 * np.cos(
                        L * np.sqrt(-k1)) - 12488 * h ** 2 * (-k1) ** (7 / 2) * (k1 + h ** 2) ** (25 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) + 1944 * h ** 2 * (-k1) ** (7 / 2) * (
                                           k1 + h ** 2) ** 13 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 1944 * h ** 2 * (-k1) ** (7 / 2) * (k1 + h ** 2) ** 13 * np.cos(
                        L * np.sqrt(-k1)) + 972 * h ** 2 * (-k1) ** (5 / 2) * (k1 + h ** 2) ** (27 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) - 92 * h ** 2 * (-k1) ** (5 / 2) * (
                                           k1 + h ** 2) ** 14 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 92 * h ** 2 * (-k1) ** (5 / 2) * (k1 + h ** 2) ** 14 * np.cos(
                        L * np.sqrt(-k1)) - 46 * h ** 2 * (-k1) ** (3 / 2) * (k1 + h ** 2) ** (29 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) + 2 * h ** 2 * (-k1) ** (3 / 2) * (
                                           k1 + h ** 2) ** 15 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) - 2 * h ** 2 * (-k1) ** (3 / 2) * (k1 + h ** 2) ** 15 * np.cos(
                        L * np.sqrt(-k1)) + h ** 2 * np.sqrt(-k1) * (k1 + h ** 2) ** (31 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) - h ** 2 * (k1 + h ** 2) ** 16 * np.sin(
                        L * np.sqrt(-k1)) - (k1 + h ** 2) ** 17 * np.sin(L * np.sqrt(-k1)) / 2) / (
                                         np.sqrt(-k1) * (k1 + h ** 2) ** 3 * (
                                         33554432 * k1 ** 14 + 142606336 * k1 ** 13 * (
                                         k1 + h ** 2) + 278921216 * k1 ** 12 * (
                                                 k1 + h ** 2) ** 2 + 332922880 * k1 ** 11 * (
                                                 k1 + h ** 2) ** 3 + 271056896 * k1 ** 10 * (
                                                 k1 + h ** 2) ** 4 + 159318016 * k1 ** 9 * (
                                                 k1 + h ** 2) ** 5 + 69746688 * k1 ** 8 * (
                                                 k1 + h ** 2) ** 6 + 23113728 * k1 ** 7 * (
                                                 k1 + h ** 2) ** 7 + 5829120 * k1 ** 6 * (
                                                 k1 + h ** 2) ** 8 + 1113728 * k1 ** 5 * (
                                                 k1 + h ** 2) ** 9 + 158752 * k1 ** 4 * (
                                                 k1 + h ** 2) ** 10 + 16376 * k1 ** 3 * (
                                                 k1 + h ** 2) ** 11 + 1156 * k1 ** 2 * (
                                                 k1 + h ** 2) ** 12 + 50 * k1 * (k1 + h ** 2) ** 13 + (
                                                 k1 + h ** 2) ** 14))
                    T[3, 0, 2] = k1 * h * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) / np.sqrt(
                        k1 + h ** 2) - h * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2))
                    T[3, 0, 3] = h * (
                            -np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) + np.sqrt(k1 + h ** 2) * (
                            np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.cos(L * np.sqrt(-k1))) / np.sqrt(k1 + h ** 2)
                    T[3, 1, 2] = -h * (
                            k1 * np.sqrt(k1 + h ** 2) * (np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.cos(L * np.sqrt(-k1)) + np.sqrt(
                        -k1) * (k1 + h ** 2) * np.sin(L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2))) / (k1 + h ** 2) ** (
                                         3 / 2)
                    T[3, 1, 3] = h * (np.sqrt(-k1) * np.sqrt(k1 + h ** 2) * (np.cos(L * np.sqrt(k1 + h ** 2)) - 1) * np.sin(
                        L * np.sqrt(-k1)) + (k1 + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1))) / (
                                         k1 + h ** 2) ** (3 / 2)
                    T[3, 2, 4] = (-k1 * (-2 * L * h ** 2 * (k1 + h ** 2) ** (3 / 2) + L * (k1 + h ** 2) ** (
                            5 / 2) + 2 * h ** 2 * (k1 + h ** 2) * np.sin(L * np.sqrt(k1 + h ** 2))) * np.cos(
                        L * np.sqrt(-k1)) + np.sqrt(-k1) * (k1 + h ** 2) ** (3 / 2) * (
                                          k1 + 2 * h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) - h ** 2) * np.sin(
                        L * np.sqrt(-k1))) / (2 * (k1 + h ** 2) ** (5 / 2))
                    T[3, 3, 4] = (-16777216 * k1 ** 15 * L * (k1 + h ** 2) ** 3 * np.sin(
                        L * np.sqrt(-k1)) - 16777216 * k1 ** 15 * h ** 2 * (k1 + h ** 2) ** 2 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 16777216 * k1 ** 15 * h ** 2 * (k1 + h ** 2) ** 2 * np.sin(
                        L * np.sqrt(-k1)) - 71303168 * k1 ** 14 * L * (k1 + h ** 2) ** 4 * np.sin(
                        L * np.sqrt(-k1)) + 8388608 * k1 ** 14 * h ** 2 * (k1 + h ** 2) ** (5 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 67108864 * k1 ** 14 * h ** 2 * (
                                          k1 + h ** 2) ** 3 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 75497472 * k1 ** 14 * h ** 2 * (k1 + h ** 2) ** 3 * np.sin(
                        L * np.sqrt(-k1)) - 139460608 * k1 ** 13 * L * (k1 + h ** 2) ** 5 * np.sin(
                        L * np.sqrt(-k1)) + 33554432 * k1 ** 13 * h ** 2 * (k1 + h ** 2) ** (7 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 122683392 * k1 ** 13 * h ** 2 * (
                                          k1 + h ** 2) ** 4 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 156237824 * k1 ** 13 * h ** 2 * (k1 + h ** 2) ** 4 * np.sin(
                        L * np.sqrt(-k1)) - 166461440 * k1 ** 12 * L * (k1 + h ** 2) ** 6 * np.sin(
                        L * np.sqrt(-k1)) + 61341696 * k1 ** 12 * h ** 2 * (k1 + h ** 2) ** (9 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 135790592 * k1 ** 12 * h ** 2 * (
                                          k1 + h ** 2) ** 5 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 197132288 * k1 ** 12 * h ** 2 * (k1 + h ** 2) ** 5 * np.sin(
                        L * np.sqrt(-k1)) - 135528448 * k1 ** 11 * L * (k1 + h ** 2) ** 7 * np.sin(
                        L * np.sqrt(-k1)) + 67895296 * k1 ** 11 * h ** 2 * (k1 + h ** 2) ** (11 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 101580800 * k1 ** 11 * h ** 2 * (
                                          k1 + h ** 2) ** 6 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 169476096 * k1 ** 11 * h ** 2 * (k1 + h ** 2) ** 6 * np.sin(
                        L * np.sqrt(-k1)) - 79659008 * k1 ** 10 * L * (k1 + h ** 2) ** 8 * np.sin(
                        L * np.sqrt(-k1)) + 50790400 * k1 ** 10 * h ** 2 * (k1 + h ** 2) ** (13 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 54263808 * k1 ** 10 * h ** 2 * (
                                          k1 + h ** 2) ** 7 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 105054208 * k1 ** 10 * h ** 2 * (k1 + h ** 2) ** 7 * np.sin(
                        L * np.sqrt(-k1)) - 34873344 * k1 ** 9 * L * (k1 + h ** 2) ** 9 * np.sin(
                        L * np.sqrt(-k1)) + 27131904 * k1 ** 9 * h ** 2 * (k1 + h ** 2) ** (15 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 21307392 * k1 ** 9 * h ** 2 * (
                                          k1 + h ** 2) ** 8 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 48439296 * k1 ** 9 * h ** 2 * (k1 + h ** 2) ** 8 * np.sin(
                        L * np.sqrt(-k1)) - 11556864 * k1 ** 8 * L * (k1 + h ** 2) ** 10 * np.sin(
                        L * np.sqrt(-k1)) + 10653696 * k1 ** 8 * h ** 2 * (k1 + h ** 2) ** (17 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 6230016 * k1 ** 8 * h ** 2 * (
                                          k1 + h ** 2) ** 9 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 16883712 * k1 ** 8 * h ** 2 * (k1 + h ** 2) ** 9 * np.sin(
                        L * np.sqrt(-k1)) - 2914560 * k1 ** 7 * L * (k1 + h ** 2) ** 11 * np.sin(
                        L * np.sqrt(-k1)) + 3115008 * k1 ** 7 * h ** 2 * (k1 + h ** 2) ** (19 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 1357056 * k1 ** 7 * h ** 2 * (
                                          k1 + h ** 2) ** 10 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 4472064 * k1 ** 7 * h ** 2 * (k1 + h ** 2) ** 10 * np.sin(
                        L * np.sqrt(-k1)) - 556864 * k1 ** 6 * L * (k1 + h ** 2) ** 12 * np.sin(
                        L * np.sqrt(-k1)) + 678528 * k1 ** 6 * h ** 2 * (k1 + h ** 2) ** (21 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 217600 * k1 ** 6 * h ** 2 * (
                                          k1 + h ** 2) ** 11 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 896128 * k1 ** 6 * h ** 2 * (k1 + h ** 2) ** 11 * np.sin(
                        L * np.sqrt(-k1)) - 79376 * k1 ** 5 * L * (k1 + h ** 2) ** 13 * np.sin(
                        L * np.sqrt(-k1)) + 108800 * k1 ** 5 * h ** 2 * (k1 + h ** 2) ** (23 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 24976 * k1 ** 5 * h ** 2 * (
                                          k1 + h ** 2) ** 12 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 133776 * k1 ** 5 * h ** 2 * (k1 + h ** 2) ** 12 * np.sin(
                        L * np.sqrt(-k1)) - 8188 * k1 ** 4 * L * (k1 + h ** 2) ** 14 * np.sin(
                        L * np.sqrt(-k1)) + 12488 * k1 ** 4 * h ** 2 * (k1 + h ** 2) ** (25 / 2) * np.sin(
                        L * np.sqrt(-k1)) * np.sin(L * np.sqrt(k1 + h ** 2)) - 1944 * k1 ** 4 * h ** 2 * (
                                          k1 + h ** 2) ** 13 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 14432 * k1 ** 4 * h ** 2 * (k1 + h ** 2) ** 13 * np.sin(
                        L * np.sqrt(-k1)) - 578 * k1 ** 3 * L * (k1 + h ** 2) ** 15 * np.sin(
                        L * np.sqrt(-k1)) + 972 * k1 ** 3 * h ** 2 * (k1 + h ** 2) ** (27 / 2) * np.sin(L * np.sqrt(-k1)) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) - 92 * k1 ** 3 * h ** 2 * (k1 + h ** 2) ** 14 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 1064 * k1 ** 3 * h ** 2 * (k1 + h ** 2) ** 14 * np.sin(
                        L * np.sqrt(-k1)) - 25 * k1 ** 2 * L * (k1 + h ** 2) ** 16 * np.sin(
                        L * np.sqrt(-k1)) + 46 * k1 ** 2 * h ** 2 * (k1 + h ** 2) ** (29 / 2) * np.sin(L * np.sqrt(-k1)) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) - 2 * k1 ** 2 * h ** 2 * (k1 + h ** 2) ** 15 * np.sin(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + 48 * k1 ** 2 * h ** 2 * (k1 + h ** 2) ** 15 * np.sin(
                        L * np.sqrt(-k1)) - k1 * L * (k1 + h ** 2) ** 17 * np.sin(L * np.sqrt(-k1)) / 2 + k1 * h ** 2 * (
                                          k1 + h ** 2) ** (31 / 2) * np.sin(L * np.sqrt(-k1)) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) + k1 * h ** 2 * (k1 + h ** 2) ** 16 * np.sin(
                        L * np.sqrt(-k1)) - 33554432 * L * h ** 2 * (-k1) ** (31 / 2) * (k1 + h ** 2) ** 2 * np.cos(
                        L * np.sqrt(-k1)) + 142606336 * L * h ** 2 * (-k1) ** (29 / 2) * (k1 + h ** 2) ** 3 * np.cos(
                        L * np.sqrt(-k1)) - 278921216 * L * h ** 2 * (-k1) ** (27 / 2) * (k1 + h ** 2) ** 4 * np.cos(
                        L * np.sqrt(-k1)) + 332922880 * L * h ** 2 * (-k1) ** (25 / 2) * (k1 + h ** 2) ** 5 * np.cos(
                        L * np.sqrt(-k1)) - 271056896 * L * h ** 2 * (-k1) ** (23 / 2) * (k1 + h ** 2) ** 6 * np.cos(
                        L * np.sqrt(-k1)) + 159318016 * L * h ** 2 * (-k1) ** (21 / 2) * (k1 + h ** 2) ** 7 * np.cos(
                        L * np.sqrt(-k1)) - 69746688 * L * h ** 2 * (-k1) ** (19 / 2) * (k1 + h ** 2) ** 8 * np.cos(
                        L * np.sqrt(-k1)) + 23113728 * L * h ** 2 * (-k1) ** (17 / 2) * (k1 + h ** 2) ** 9 * np.cos(
                        L * np.sqrt(-k1)) - 5829120 * L * h ** 2 * (-k1) ** (15 / 2) * (k1 + h ** 2) ** 10 * np.cos(
                        L * np.sqrt(-k1)) + 1113728 * L * h ** 2 * (-k1) ** (13 / 2) * (k1 + h ** 2) ** 11 * np.cos(
                        L * np.sqrt(-k1)) - 158752 * L * h ** 2 * (-k1) ** (11 / 2) * (k1 + h ** 2) ** 12 * np.cos(
                        L * np.sqrt(-k1)) + 16376 * L * h ** 2 * (-k1) ** (9 / 2) * (k1 + h ** 2) ** 13 * np.cos(
                        L * np.sqrt(-k1)) - 1156 * L * h ** 2 * (-k1) ** (7 / 2) * (k1 + h ** 2) ** 14 * np.cos(
                        L * np.sqrt(-k1)) + 50 * L * h ** 2 * (-k1) ** (5 / 2) * (k1 + h ** 2) ** 15 * np.cos(
                        L * np.sqrt(-k1)) - L * h ** 2 * (-k1) ** (3 / 2) * (k1 + h ** 2) ** 16 * np.cos(
                        L * np.sqrt(-k1)) + 33554432 * h ** 2 * (-k1) ** (31 / 2) * (k1 + h ** 2) ** (3 / 2) * np.sin(
                        L * np.sqrt(k1 + h ** 2)) * np.cos(L * np.sqrt(-k1)) - 150994944 * h ** 2 * (-k1) ** (29 / 2) * (
                                          k1 + h ** 2) ** (5 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) - 16777216 * h ** 2 * (-k1) ** (29 / 2) * (k1 + h ** 2) ** 2 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) + 16777216 * h ** 2 * (-k1) ** (29 / 2) * (
                                          k1 + h ** 2) ** 2 * np.cos(L * np.sqrt(-k1)) + 312475648 * h ** 2 * (-k1) ** (
                                          27 / 2) * (k1 + h ** 2) ** (7 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) + 75497472 * h ** 2 * (-k1) ** (27 / 2) * (k1 + h ** 2) ** 3 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) - 75497472 * h ** 2 * (-k1) ** (27 / 2) * (
                                          k1 + h ** 2) ** 3 * np.cos(L * np.sqrt(-k1)) - 394264576 * h ** 2 * (-k1) ** (
                                          25 / 2) * (k1 + h ** 2) ** (9 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) - 156237824 * h ** 2 * (-k1) ** (25 / 2) * (k1 + h ** 2) ** 4 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) + 156237824 * h ** 2 * (-k1) ** (25 / 2) * (
                                          k1 + h ** 2) ** 4 * np.cos(L * np.sqrt(-k1)) + 338952192 * h ** 2 * (-k1) ** (
                                          23 / 2) * (k1 + h ** 2) ** (11 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) + 197132288 * h ** 2 * (-k1) ** (23 / 2) * (k1 + h ** 2) ** 5 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) - 197132288 * h ** 2 * (-k1) ** (23 / 2) * (
                                          k1 + h ** 2) ** 5 * np.cos(L * np.sqrt(-k1)) - 210108416 * h ** 2 * (-k1) ** (
                                          21 / 2) * (k1 + h ** 2) ** (13 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) - 169476096 * h ** 2 * (-k1) ** (21 / 2) * (k1 + h ** 2) ** 6 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) + 169476096 * h ** 2 * (-k1) ** (21 / 2) * (
                                          k1 + h ** 2) ** 6 * np.cos(L * np.sqrt(-k1)) + 96878592 * h ** 2 * (-k1) ** (
                                          19 / 2) * (k1 + h ** 2) ** (15 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) + 105054208 * h ** 2 * (-k1) ** (19 / 2) * (k1 + h ** 2) ** 7 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) - 105054208 * h ** 2 * (-k1) ** (19 / 2) * (
                                          k1 + h ** 2) ** 7 * np.cos(L * np.sqrt(-k1)) - 33767424 * h ** 2 * (-k1) ** (
                                          17 / 2) * (k1 + h ** 2) ** (17 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) - 48439296 * h ** 2 * (-k1) ** (17 / 2) * (k1 + h ** 2) ** 8 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) + 48439296 * h ** 2 * (-k1) ** (17 / 2) * (
                                          k1 + h ** 2) ** 8 * np.cos(L * np.sqrt(-k1)) + 8944128 * h ** 2 * (-k1) ** (
                                          15 / 2) * (k1 + h ** 2) ** (19 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) + 16883712 * h ** 2 * (-k1) ** (15 / 2) * (k1 + h ** 2) ** 9 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) - 16883712 * h ** 2 * (-k1) ** (15 / 2) * (
                                          k1 + h ** 2) ** 9 * np.cos(L * np.sqrt(-k1)) - 1792256 * h ** 2 * (-k1) ** (
                                          13 / 2) * (k1 + h ** 2) ** (21 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) - 4472064 * h ** 2 * (-k1) ** (13 / 2) * (k1 + h ** 2) ** 10 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) + 4472064 * h ** 2 * (-k1) ** (13 / 2) * (
                                          k1 + h ** 2) ** 10 * np.cos(L * np.sqrt(-k1)) + 267552 * h ** 2 * (-k1) ** (
                                          11 / 2) * (k1 + h ** 2) ** (23 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) + 896128 * h ** 2 * (-k1) ** (11 / 2) * (k1 + h ** 2) ** 11 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) - 896128 * h ** 2 * (-k1) ** (11 / 2) * (
                                          k1 + h ** 2) ** 11 * np.cos(L * np.sqrt(-k1)) - 28864 * h ** 2 * (-k1) ** (
                                          9 / 2) * (k1 + h ** 2) ** (25 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) - 133776 * h ** 2 * (-k1) ** (9 / 2) * (k1 + h ** 2) ** 12 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) + 133776 * h ** 2 * (-k1) ** (9 / 2) * (
                                          k1 + h ** 2) ** 12 * np.cos(L * np.sqrt(-k1)) + 2128 * h ** 2 * (-k1) ** (
                                          7 / 2) * (k1 + h ** 2) ** (27 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) + 14432 * h ** 2 * (-k1) ** (7 / 2) * (k1 + h ** 2) ** 13 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) - 14432 * h ** 2 * (-k1) ** (7 / 2) * (
                                          k1 + h ** 2) ** 13 * np.cos(L * np.sqrt(-k1)) - 96 * h ** 2 * (-k1) ** (
                                          5 / 2) * (k1 + h ** 2) ** (29 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) - 1064 * h ** 2 * (-k1) ** (5 / 2) * (k1 + h ** 2) ** 14 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) + 1064 * h ** 2 * (-k1) ** (5 / 2) * (
                                          k1 + h ** 2) ** 14 * np.cos(L * np.sqrt(-k1)) + 2 * h ** 2 * (-k1) ** (
                                          3 / 2) * (k1 + h ** 2) ** (31 / 2) * np.sin(L * np.sqrt(k1 + h ** 2)) * np.cos(
                        L * np.sqrt(-k1)) + 48 * h ** 2 * (-k1) ** (3 / 2) * (k1 + h ** 2) ** 15 * np.cos(
                        L * np.sqrt(-k1)) * np.cos(L * np.sqrt(k1 + h ** 2)) - 48 * h ** 2 * (-k1) ** (3 / 2) * (
                                          k1 + h ** 2) ** 15 * np.cos(L * np.sqrt(-k1)) - h ** 2 * np.sqrt(-k1) * (
                                          k1 + h ** 2) ** 16 * np.cos(L * np.sqrt(-k1)) * np.cos(
                        L * np.sqrt(k1 + h ** 2)) + h ** 2 * np.sqrt(-k1) * (k1 + h ** 2) ** 16 * np.cos(L * np.sqrt(-k1))) / (
                                         np.sqrt(-k1) * (k1 + h ** 2) ** 3 * (
                                         33554432 * k1 ** 14 + 142606336 * k1 ** 13 * (
                                         k1 + h ** 2) + 278921216 * k1 ** 12 * (
                                                 k1 + h ** 2) ** 2 + 332922880 * k1 ** 11 * (
                                                 k1 + h ** 2) ** 3 + 271056896 * k1 ** 10 * (
                                                 k1 + h ** 2) ** 4 + 159318016 * k1 ** 9 * (
                                                 k1 + h ** 2) ** 5 + 69746688 * k1 ** 8 * (
                                                 k1 + h ** 2) ** 6 + 23113728 * k1 ** 7 * (
                                                 k1 + h ** 2) ** 7 + 5829120 * k1 ** 6 * (
                                                 k1 + h ** 2) ** 8 + 1113728 * k1 ** 5 * (
                                                 k1 + h ** 2) ** 9 + 158752 * k1 ** 4 * (
                                                 k1 + h ** 2) ** 10 + 16376 * k1 ** 3 * (
                                                 k1 + h ** 2) ** 11 + 1156 * k1 ** 2 * (
                                                 k1 + h ** 2) ** 12 + 50 * k1 * (k1 + h ** 2) ** 13 + (
                                                 k1 + h ** 2) ** 14))
            if h ** 2 + k1 < 0:
                T[0, 0, 0] = h * (4 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) ** 2 + 4 * k1 * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - 7 * k1 + 2 * h ** 2 * np.cos(L * np.sqrt(k1 + h ** 2)) ** 2 + 2 * h ** 2 * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - 3 * h ** 2 + (k1 + h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2)) ** 2 - 2 * (
                                          k1 + h ** 2) * np.cosh(L * np.sqrt(-k1 - h ** 2))) / (6 * (k1 + h ** 2))
                T[0, 0, 1] = -h * (4 * k1 * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 5 * k1 + 2 * h ** 2 * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - 3 * h ** 2 + (k1 + h ** 2) * np.cosh(L * np.sqrt(-k1 - h ** 2))) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) / (3 * (-k1 - h ** 2) ** (3 / 2))
                T[0, 0, 4] = (-6 * L * h ** 2 * (-k1 - h ** 2) ** (7 / 2) * (2 * k1 + h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - 3 * L * (-k1 - h ** 2) ** (9 / 2) * (k1 + 2 * h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 2 * h ** 2 * (k1 + h ** 2) ** 4 * (
                                      -np.sinh(L * np.sqrt(-k1 - h ** 2)) ** 2 + 2 * np.cosh(
                                  L * np.sqrt(-k1 - h ** 2)) - 2) + 4 * h ** 2 * (k1 + h ** 2) ** 3 * (
                                      -2 * k1 * np.sinh(L * np.sqrt(-k1 - h ** 2)) ** 2 - 2 * k1 * np.cosh(
                                  L * np.sqrt(-k1 - h ** 2)) + 2 * k1 - h ** 2 * np.sinh(
                                  L * np.sqrt(-k1 - h ** 2)) ** 2 - h ** 2 * np.cosh(
                                  L * np.sqrt(-k1 - h ** 2)) + h ** 2)) / (6 * (k1 + h ** 2) ** 5)
                T[0, 1, 1] = h * (-4 * k1 * np.sinh(L * np.sqrt(-k1 - h ** 2)) ** 2 + 8 * k1 * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - 7 * k1 - 2 * h ** 2 * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) ** 2 + 4 * h ** 2 * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 3 * h ** 2 - (
                                          k1 + h ** 2) * np.sinh(L * np.sqrt(-k1 - h ** 2)) ** 2 - (k1 + h ** 2) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2))) / (6 * (k1 + h ** 2) ** 2)
                T[0, 1, 4] = (-L * h ** 2 * (k1 + h ** 2) ** 5 * (2 * k1 + h ** 2) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - L * (k1 + h ** 2) ** 6 * (k1 + 2 * h ** 2) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) / 2 - h ** 2 * (-k1 - h ** 2) ** (9 / 2) * (
                                      4 * k1 * np.cosh(L * np.sqrt(-k1 - h ** 2)) + 2 * k1 + 2 * h ** 2 * np.cosh(
                                  L * np.sqrt(-k1 - h ** 2)) + h ** 2) * np.sinh(L * np.sqrt(-k1 - h ** 2)) / 3 + (
                                      -k1 - h ** 2) ** (11 / 2) * (
                                      3 * k1 - 2 * h ** 2 * np.cosh(L * np.sqrt(-k1 - h ** 2)) + 8 * h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) / 6) / (k1 + h ** 2) ** 7
                T[0, 2, 2] = k1 * h * (-np.cos(2 * L * np.sqrt(-k1)) + np.cosh(L * np.sqrt(-k1 - h ** 2))) / (2 * (5 * k1 + h ** 2))
                T[0, 3, 3] = h * (np.cosh(L * np.sqrt(-k1 - h ** 2)) - 1) / (2 * (k1 + h ** 2))
                T[0, 4, 4] = h * (-6 * L * h ** 2 * (-k1 - h ** 2) ** (15 / 2) * (2 * k1 + h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 3 * L * (-k1 - h ** 2) ** (17 / 2) * (k1 + 2 * h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 2 * h ** 2 * (k1 + h ** 2) ** 7 * (
                                          6 * k1 * np.cos(L * np.sqrt(k1 + h ** 2)) ** 6 - 6 * k1 * np.sinh(
                                      L * np.sqrt(-k1 - h ** 2)) ** 6 - 18 * k1 * np.sinh(
                                      L * np.sqrt(-k1 - h ** 2)) ** 4 - 20 * k1 * np.sinh(
                                      L * np.sqrt(-k1 - h ** 2)) ** 2 - 8 * k1 * np.cosh(
                                      L * np.sqrt(-k1 - h ** 2)) + 2 * k1 + 3 * h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) ** 6 - 3 * h ** 2 * np.sinh(
                                      L * np.sqrt(-k1 - h ** 2)) ** 6 - 9 * h ** 2 * np.sinh(
                                      L * np.sqrt(-k1 - h ** 2)) ** 4 - 10 * h ** 2 * np.sinh(
                                      L * np.sqrt(-k1 - h ** 2)) ** 2 - 4 * h ** 2 * np.cosh(
                                      L * np.sqrt(-k1 - h ** 2)) + h ** 2) + 6 * (k1 + h ** 2) ** 9 * (
                                          np.cosh(L * np.sqrt(-k1 - h ** 2)) - 1) + (k1 + h ** 2) ** 8 * (
                                          -6 * k1 * np.cosh(L * np.sqrt(-k1 - h ** 2)) + 6 * k1 + h ** 2 * np.cos(
                                      L * np.sqrt(k1 + h ** 2)) ** 2 - 14 * h ** 2 * np.cosh(
                                      L * np.sqrt(-k1 - h ** 2)) + 13 * h ** 2)) / (6 * (k1 + h ** 2) ** 10)
                T[1, 0, 0] = -h * (-2 * k1 * np.sqrt(-k1 - h ** 2) * np.sinh(L * np.sqrt(-k1 - h ** 2)) - 2 * k1 * np.sqrt(
                    -k1 - h ** 2) * np.sinh(2 * L * np.sqrt(-k1 - h ** 2)) - h ** 2 * np.sqrt(-k1 - h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - h ** 2 * np.sqrt(-k1 - h ** 2) * np.sinh(2 * L * np.sqrt(-k1 - h ** 2)) - (
                                           -k1 - h ** 2) ** (3 / 2) * np.sinh(L * np.sqrt(-k1 - h ** 2)) + (
                                           -k1 - h ** 2) ** (3 / 2) * np.sinh(2 * L * np.sqrt(-k1 - h ** 2)) / 2) / (
                                     3 * (k1 + h ** 2))
                T[1, 0, 1] = h * ((5 * k1 + 3 * h ** 2) * np.sinh(L * np.sqrt(-k1 - h ** 2)) ** 2 + (
                        4 * k1 * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 5 * k1 + 2 * h ** 2 * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - 3 * h ** 2 + (k1 + h ** 2) * np.cosh(L * np.sqrt(-k1 - h ** 2))) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2))) / (3 * (k1 + h ** 2))
                T[1, 0, 4] = -(
                        L * h ** 2 * (k1 + h ** 2) ** 4 * (2 * k1 + h ** 2) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - L * (
                        k1 + h ** 2) ** 5 * (k1 + 2 * h ** 2) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) / 2 + 2 * h ** 2 * (-k1 - h ** 2) ** (9 / 2) * (
                                np.cosh(L * np.sqrt(-k1 - h ** 2)) - 1) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) / 3 + h ** 2 * (-k1 - h ** 2) ** (7 / 2) * (2 * k1 + h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - 2 * h ** 2 * (-k1 - h ** 2) ** (7 / 2) * (
                                4 * k1 * np.cosh(L * np.sqrt(-k1 - h ** 2)) + 2 * k1 + 2 * h ** 2 * np.cosh(
                            L * np.sqrt(-k1 - h ** 2)) + h ** 2) * np.sinh(L * np.sqrt(-k1 - h ** 2)) / 3 + (
                                -k1 - h ** 2) ** (9 / 2) * (k1 + 2 * h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) / 2) / (k1 + h ** 2) ** 5
                T[1, 1, 1] = -h * (8 * k1 * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 7 * k1 + 4 * h ** 2 * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - 3 * h ** 2 + 2 * (k1 + h ** 2) * np.cosh(L * np.sqrt(-k1 - h ** 2))) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) / (6 * (-k1 - h ** 2) ** (3 / 2))
                T[1, 1, 4] = (12 * k1 * L * h ** 2 * (-k1 - h ** 2) ** (11 / 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - 3 * k1 * L * (-k1 - h ** 2) ** (13 / 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 16 * k1 * h ** 2 * (k1 + h ** 2) ** 5 * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) ** 2 - 8 * k1 * h ** 2 * (k1 + h ** 2) ** 5 * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) + 8 * k1 * h ** 2 * (k1 + h ** 2) ** 5 + 6 * L * h ** 4 * (
                                      -k1 - h ** 2) ** (11 / 2) * np.sinh(L * np.sqrt(-k1 - h ** 2)) - 6 * L * h ** 2 * (
                                      -k1 - h ** 2) ** (13 / 2) * np.sinh(L * np.sqrt(-k1 - h ** 2)) + 8 * h ** 4 * (
                                      k1 + h ** 2) ** 5 * np.sinh(L * np.sqrt(-k1 - h ** 2)) ** 2 - 4 * h ** 4 * (
                                      k1 + h ** 2) ** 5 * np.cosh(L * np.sqrt(-k1 - h ** 2)) + 4 * h ** 4 * (
                                      k1 + h ** 2) ** 5 - 4 * h ** 2 * (k1 + h ** 2) ** 6 * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) ** 2 + 2 * h ** 2 * (k1 + h ** 2) ** 6 * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - 2 * h ** 2 * (k1 + h ** 2) ** 6) / (6 * (k1 + h ** 2) ** 7)
                T[1, 2, 2] = k1 * h * (2 * np.sqrt(-k1) * np.sin(2 * L * np.sqrt(-k1)) + np.sqrt(-k1 - h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2))) / (2 * (5 * k1 + h ** 2))
                T[1, 3, 3] = -h * np.sinh(L * np.sqrt(-k1 - h ** 2)) / (2 * np.sqrt(-k1 - h ** 2))
                T[1, 4, 4] = -h * (6 * L * h ** 2 * (k1 + h ** 2) ** 8 * (2 * k1 + h ** 2) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) + 3 * L * (k1 + h ** 2) ** 9 * (k1 + 2 * h ** 2) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) + 6 * h ** 2 * (-k1 - h ** 2) ** (15 / 2) * (2 * k1 + h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + h ** 2 * (k1 + h ** 2) ** 7 * (
                                           -72 * k1 * np.sqrt(-k1 - h ** 2) * np.cos(L * np.sqrt(k1 + h ** 2)) ** 5 * np.sinh(
                                       L * np.sqrt(-k1 - h ** 2)) + 72 * k1 * np.sqrt(-k1 - h ** 2) * np.sinh(
                                       L * np.sqrt(-k1 - h ** 2)) ** 5 * np.cosh(
                                       L * np.sqrt(-k1 - h ** 2)) + 144 * k1 * np.sqrt(-k1 - h ** 2) * np.sinh(
                                       L * np.sqrt(-k1 - h ** 2)) ** 3 * np.cosh(L * np.sqrt(-k1 - h ** 2)) + 80 * k1 * np.sqrt(
                                       -k1 - h ** 2) * np.sinh(L * np.sqrt(-k1 - h ** 2)) * np.cosh(
                                       L * np.sqrt(-k1 - h ** 2)) + 16 * k1 * np.sqrt(-k1 - h ** 2) * np.sinh(
                                       L * np.sqrt(-k1 - h ** 2)) - 36 * h ** 2 * np.sqrt(-k1 - h ** 2) * np.cos(
                                       L * np.sqrt(k1 + h ** 2)) ** 5 * np.sinh(
                                       L * np.sqrt(-k1 - h ** 2)) + 36 * h ** 2 * np.sqrt(-k1 - h ** 2) * np.sinh(
                                       L * np.sqrt(-k1 - h ** 2)) ** 5 * np.cosh(
                                       L * np.sqrt(-k1 - h ** 2)) + 72 * h ** 2 * np.sqrt(-k1 - h ** 2) * np.sinh(
                                       L * np.sqrt(-k1 - h ** 2)) ** 3 * np.cosh(
                                       L * np.sqrt(-k1 - h ** 2)) + 40 * h ** 2 * np.sqrt(-k1 - h ** 2) * np.sinh(
                                       L * np.sqrt(-k1 - h ** 2)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) + 8 * h ** 2 * np.sqrt(
                                       -k1 - h ** 2) * np.sinh(L * np.sqrt(-k1 - h ** 2))) + 6 * (-k1 - h ** 2) ** (
                                           19 / 2) * np.sinh(L * np.sqrt(-k1 - h ** 2)) - 3 * (-k1 - h ** 2) ** (
                                           17 / 2) * (k1 + 2 * h ** 2) * np.sinh(L * np.sqrt(-k1 - h ** 2)) + (
                                           k1 + h ** 2) ** 8 * (6 * k1 * np.sqrt(-k1 - h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 14 * h ** 2 * np.sqrt(-k1 - h ** 2) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - h ** 2 * np.sqrt(-k1 - h ** 2) * np.sinh(2 * L * np.sqrt(-k1 - h ** 2)))) / (
                                     6 * (k1 + h ** 2) ** 10)
                T[2, 0, 2] = -h * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1 - h ** 2)) / np.sqrt(-k1 - h ** 2)
                T[2, 0, 3] = -h * (
                        -np.sqrt(-k1) * np.cos(L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1 - h ** 2)) + np.sqrt(-k1 - h ** 2) * np.sin(
                    L * np.sqrt(-k1))) / (np.sqrt(-k1) * np.sqrt(-k1 - h ** 2))
                T[2, 1, 2] = h * np.sqrt(-k1) * (np.cosh(L * np.sqrt(-k1 - h ** 2)) - 1) * np.sin(L * np.sqrt(-k1)) / (k1 + h ** 2)
                T[2, 1, 3] = -h * (np.cosh(L * np.sqrt(-k1 - h ** 2)) - 1) * np.cos(L * np.sqrt(-k1)) / (k1 + h ** 2)
                T[2, 2, 4] = np.sqrt(-k1) * (
                        2 * L * h ** 2 * (-k1 - h ** 2) ** (3 / 2) + L * (-k1 - h ** 2) ** (5 / 2) + 2 * h ** 2 * (
                        k1 + h ** 2) * np.sinh(L * np.sqrt(-k1 - h ** 2))) * np.sin(L * np.sqrt(-k1)) / (
                                     2 * (-k1 - h ** 2) ** (5 / 2))
                T[2, 3, 4] = (-33554432 * k1 ** 15 * h ** 2 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - 125829120 * k1 ** 14 * L * h ** 2 * (-k1 - h ** 2) ** (3 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 117440512 * k1 ** 14 * h ** 2 * (k1 + h ** 2) * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - 16777216 * k1 ** 14 * (-k1 - h ** 2) ** (3 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 216006656 * k1 ** 13 * L * h ** 2 * (-k1 - h ** 2) ** (5 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 58720256 * k1 ** 13 * h ** 2 * (-k1 - h ** 2) ** (3 / 2) * np.sin(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 67108864 * k1 ** 13 * h ** 2 * (-k1 - h ** 2) ** (
                                      3 / 2) * np.sin(L * np.sqrt(-k1)) - 186646528 * k1 ** 13 * h ** 2 * (
                                      k1 + h ** 2) ** 2 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 62914560 * k1 ** 13 * (-k1 - h ** 2) ** (5 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 224919552 * k1 ** 12 * L * h ** 2 * (-k1 - h ** 2) ** (7 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 93323264 * k1 ** 12 * h ** 2 * (-k1 - h ** 2) ** (5 / 2) * np.sin(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) + 122683392 * k1 ** 12 * h ** 2 * (-k1 - h ** 2) ** (
                                      5 / 2) * np.sin(L * np.sqrt(-k1)) - 178257920 * k1 ** 12 * h ** 2 * (
                                      k1 + h ** 2) ** 3 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - 108003328 * k1 ** 12 * (-k1 - h ** 2) ** (7 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 158597120 * k1 ** 11 * L * h ** 2 * (-k1 - h ** 2) ** (9 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 89128960 * k1 ** 11 * h ** 2 * (-k1 - h ** 2) ** (7 / 2) * np.sin(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 135790592 * k1 ** 11 * h ** 2 * (-k1 - h ** 2) ** (
                                      7 / 2) * np.sin(L * np.sqrt(-k1)) - 114032640 * k1 ** 11 * h ** 2 * (
                                      k1 + h ** 2) ** 4 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 112459776 * k1 ** 11 * (-k1 - h ** 2) ** (9 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 80019456 * k1 ** 10 * L * h ** 2 * (-k1 - h ** 2) ** (11 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 57016320 * k1 ** 10 * h ** 2 * (-k1 - h ** 2) ** (9 / 2) * np.sin(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) + 101580800 * k1 ** 10 * h ** 2 * (-k1 - h ** 2) ** (
                                      9 / 2) * np.sin(L * np.sqrt(-k1)) - 51511296 * k1 ** 10 * h ** 2 * (
                                      k1 + h ** 2) ** 5 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - 79298560 * k1 ** 10 * (-k1 - h ** 2) ** (11 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 29736960 * k1 ** 9 * L * h ** 2 * (-k1 - h ** 2) ** (13 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 25755648 * k1 ** 9 * h ** 2 * (-k1 - h ** 2) ** (11 / 2) * np.sin(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 54263808 * k1 ** 9 * h ** 2 * (-k1 - h ** 2) ** (
                                      11 / 2) * np.sin(L * np.sqrt(-k1)) - 16859136 * k1 ** 9 * h ** 2 * (
                                      k1 + h ** 2) ** 6 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 40009728 * k1 ** 9 * (-k1 - h ** 2) ** (13 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 8245248 * k1 ** 8 * L * h ** 2 * (-k1 - h ** 2) ** (15 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 8429568 * k1 ** 8 * h ** 2 * (-k1 - h ** 2) ** (13 / 2) * np.sin(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) + 21307392 * k1 ** 8 * h ** 2 * (-k1 - h ** 2) ** (
                                      13 / 2) * np.sin(L * np.sqrt(-k1)) - 4030464 * k1 ** 8 * h ** 2 * (
                                      k1 + h ** 2) ** 7 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - 14868480 * k1 ** 8 * (-k1 - h ** 2) ** (15 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 1706496 * k1 ** 7 * L * h ** 2 * (-k1 - h ** 2) ** (17 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 2015232 * k1 ** 7 * h ** 2 * (-k1 - h ** 2) ** (15 / 2) * np.sin(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 6230016 * k1 ** 7 * h ** 2 * (-k1 - h ** 2) ** (
                                      15 / 2) * np.sin(L * np.sqrt(-k1)) - 698880 * k1 ** 7 * h ** 2 * (
                                      k1 + h ** 2) ** 8 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 4122624 * k1 ** 7 * (-k1 - h ** 2) ** (17 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 260480 * k1 ** 6 * L * h ** 2 * (-k1 - h ** 2) ** (19 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 349440 * k1 ** 6 * h ** 2 * (-k1 - h ** 2) ** (17 / 2) * np.sin(L * np.sqrt(-k1)) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) + 1357056 * k1 ** 6 * h ** 2 * (-k1 - h ** 2) ** (17 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 85760 * k1 ** 6 * h ** 2 * (k1 + h ** 2) ** 9 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - 853248 * k1 ** 6 * (-k1 - h ** 2) ** (19 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 28512 * k1 ** 5 * L * h ** 2 * (-k1 - h ** 2) ** (21 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 42880 * k1 ** 5 * h ** 2 * (-k1 - h ** 2) ** (19 / 2) * np.sin(L * np.sqrt(-k1)) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - 217600 * k1 ** 5 * h ** 2 * (-k1 - h ** 2) ** (19 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 7072 * k1 ** 5 * h ** 2 * (k1 + h ** 2) ** 10 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 130240 * k1 ** 5 * (-k1 - h ** 2) ** (21 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 2120 * k1 ** 4 * L * h ** 2 * (-k1 - h ** 2) ** (23 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 3536 * k1 ** 4 * h ** 2 * (-k1 - h ** 2) ** (21 / 2) * np.sin(L * np.sqrt(-k1)) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) + 24976 * k1 ** 4 * h ** 2 * (-k1 - h ** 2) ** (21 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 352 * k1 ** 4 * h ** 2 * (k1 + h ** 2) ** 11 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - 14256 * k1 ** 4 * (-k1 - h ** 2) ** (23 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 96 * k1 ** 3 * L * h ** 2 * (-k1 - h ** 2) ** (25 / 2) * np.sin(
                    L * np.sqrt(-k1)) + 176 * k1 ** 3 * h ** 2 * (-k1 - h ** 2) ** (23 / 2) * np.sin(L * np.sqrt(-k1)) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - 1944 * k1 ** 3 * h ** 2 * (-k1 - h ** 2) ** (23 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 8 * k1 ** 3 * h ** 2 * (k1 + h ** 2) ** 12 * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 1060 * k1 ** 3 * (-k1 - h ** 2) ** (25 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 2 * k1 ** 2 * L * h ** 2 * (-k1 - h ** 2) ** (27 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 4 * k1 ** 2 * h ** 2 * (-k1 - h ** 2) ** (25 / 2) * np.sin(L * np.sqrt(-k1)) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) + 92 * k1 ** 2 * h ** 2 * (-k1 - h ** 2) ** (25 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 48 * k1 ** 2 * (-k1 - h ** 2) ** (27 / 2) * np.sin(
                    L * np.sqrt(-k1)) - 2 * k1 * h ** 2 * (-k1 - h ** 2) ** (27 / 2) * np.sin(L * np.sqrt(-k1)) + k1 * (
                                      -k1 - h ** 2) ** (29 / 2) * np.sin(L * np.sqrt(-k1)) + 16777216 * L * (-k1) ** (
                                      29 / 2) * (-k1 - h ** 2) ** (3 / 2) * np.cos(L * np.sqrt(-k1)) + 62914560 * L * (
                                  -k1) ** (27 / 2) * (-k1 - h ** 2) ** (5 / 2) * np.cos(L * np.sqrt(-k1)) + 108003328 * L * (
                                  -k1) ** (25 / 2) * (-k1 - h ** 2) ** (7 / 2) * np.cos(L * np.sqrt(-k1)) + 112459776 * L * (
                                  -k1) ** (23 / 2) * (-k1 - h ** 2) ** (9 / 2) * np.cos(L * np.sqrt(-k1)) + 79298560 * L * (
                                  -k1) ** (21 / 2) * (-k1 - h ** 2) ** (11 / 2) * np.cos(L * np.sqrt(-k1)) + 40009728 * L * (
                                  -k1) ** (19 / 2) * (-k1 - h ** 2) ** (13 / 2) * np.cos(L * np.sqrt(-k1)) + 14868480 * L * (
                                  -k1) ** (17 / 2) * (-k1 - h ** 2) ** (15 / 2) * np.cos(L * np.sqrt(-k1)) + 4122624 * L * (
                                  -k1) ** (15 / 2) * (-k1 - h ** 2) ** (17 / 2) * np.cos(L * np.sqrt(-k1)) + 853248 * L * (
                                  -k1) ** (13 / 2) * (-k1 - h ** 2) ** (19 / 2) * np.cos(L * np.sqrt(-k1)) + 130240 * L * (
                                  -k1) ** (11 / 2) * (-k1 - h ** 2) ** (21 / 2) * np.cos(L * np.sqrt(-k1)) + 14256 * L * (
                                  -k1) ** (9 / 2) * (-k1 - h ** 2) ** (23 / 2) * np.cos(L * np.sqrt(-k1)) + 1060 * L * (
                                  -k1) ** (7 / 2) * (-k1 - h ** 2) ** (25 / 2) * np.cos(L * np.sqrt(-k1)) + 48 * L * (
                                  -k1) ** (5 / 2) * (-k1 - h ** 2) ** (27 / 2) * np.cos(L * np.sqrt(-k1)) + L * (-k1) ** (
                                      3 / 2) * (-k1 - h ** 2) ** (29 / 2) * np.cos(
                    L * np.sqrt(-k1)) + 58720256 * h ** 2 * (-k1) ** (27 / 2) * (-k1 - h ** 2) ** (3 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 58720256 * h ** 2 * (-k1) ** (27 / 2) * (
                                      -k1 - h ** 2) ** (3 / 2) * np.cos(L * np.sqrt(-k1)) + 8388608 * h ** 2 * (-k1) ** (
                                      27 / 2) * (k1 + h ** 2) * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 93323264 * h ** 2 * (-k1) ** (25 / 2) * (-k1 - h ** 2) ** (5 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 93323264 * h ** 2 * (-k1) ** (25 / 2) * (
                                      -k1 - h ** 2) ** (5 / 2) * np.cos(L * np.sqrt(-k1)) - 29360128 * h ** 2 * (-k1) ** (
                                      25 / 2) * (k1 + h ** 2) ** 2 * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 89128960 * h ** 2 * (-k1) ** (23 / 2) * (-k1 - h ** 2) ** (7 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 89128960 * h ** 2 * (-k1) ** (23 / 2) * (
                                      -k1 - h ** 2) ** (7 / 2) * np.cos(L * np.sqrt(-k1)) + 46661632 * h ** 2 * (-k1) ** (
                                      23 / 2) * (k1 + h ** 2) ** 3 * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 57016320 * h ** 2 * (-k1) ** (21 / 2) * (-k1 - h ** 2) ** (9 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 57016320 * h ** 2 * (-k1) ** (21 / 2) * (
                                      -k1 - h ** 2) ** (9 / 2) * np.cos(L * np.sqrt(-k1)) - 44564480 * h ** 2 * (-k1) ** (
                                      21 / 2) * (k1 + h ** 2) ** 4 * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 25755648 * h ** 2 * (-k1) ** (19 / 2) * (-k1 - h ** 2) ** (11 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 25755648 * h ** 2 * (-k1) ** (19 / 2) * (
                                      -k1 - h ** 2) ** (11 / 2) * np.cos(L * np.sqrt(-k1)) + 28508160 * h ** 2 * (
                                  -k1) ** (19 / 2) * (k1 + h ** 2) ** 5 * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 8429568 * h ** 2 * (-k1) ** (17 / 2) * (-k1 - h ** 2) ** (13 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 8429568 * h ** 2 * (-k1) ** (17 / 2) * (
                                      -k1 - h ** 2) ** (13 / 2) * np.cos(L * np.sqrt(-k1)) - 12877824 * h ** 2 * (
                                  -k1) ** (17 / 2) * (k1 + h ** 2) ** 6 * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 2015232 * h ** 2 * (-k1) ** (15 / 2) * (-k1 - h ** 2) ** (15 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 2015232 * h ** 2 * (-k1) ** (15 / 2) * (
                                      -k1 - h ** 2) ** (15 / 2) * np.cos(L * np.sqrt(-k1)) + 4214784 * h ** 2 * (-k1) ** (
                                      15 / 2) * (k1 + h ** 2) ** 7 * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 349440 * h ** 2 * (-k1) ** (13 / 2) * (-k1 - h ** 2) ** (17 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 349440 * h ** 2 * (-k1) ** (13 / 2) * (
                                      -k1 - h ** 2) ** (17 / 2) * np.cos(L * np.sqrt(-k1)) - 1007616 * h ** 2 * (-k1) ** (
                                      13 / 2) * (k1 + h ** 2) ** 8 * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 42880 * h ** 2 * (-k1) ** (11 / 2) * (-k1 - h ** 2) ** (19 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 42880 * h ** 2 * (-k1) ** (11 / 2) * (
                                      -k1 - h ** 2) ** (19 / 2) * np.cos(L * np.sqrt(-k1)) + 174720 * h ** 2 * (-k1) ** (
                                      11 / 2) * (k1 + h ** 2) ** 9 * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 3536 * h ** 2 * (-k1) ** (9 / 2) * (-k1 - h ** 2) ** (21 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 3536 * h ** 2 * (-k1) ** (9 / 2) * (
                                      -k1 - h ** 2) ** (21 / 2) * np.cos(L * np.sqrt(-k1)) - 21440 * h ** 2 * (-k1) ** (
                                      9 / 2) * (k1 + h ** 2) ** 10 * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 176 * h ** 2 * (-k1) ** (7 / 2) * (-k1 - h ** 2) ** (23 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 176 * h ** 2 * (-k1) ** (7 / 2) * (
                                      -k1 - h ** 2) ** (23 / 2) * np.cos(L * np.sqrt(-k1)) + 1768 * h ** 2 * (-k1) ** (
                                      7 / 2) * (k1 + h ** 2) ** 11 * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 4 * h ** 2 * (-k1) ** (5 / 2) * (-k1 - h ** 2) ** (25 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - 4 * h ** 2 * (-k1) ** (5 / 2) * (-k1 - h ** 2) ** (
                                      25 / 2) * np.cos(L * np.sqrt(-k1)) - 88 * h ** 2 * (-k1) ** (5 / 2) * (
                                      k1 + h ** 2) ** 12 * np.cos(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) + 2 * h ** 2 * (-k1) ** (3 / 2) * (k1 + h ** 2) ** 13 * np.cos(
                    L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1 - h ** 2)) + 16777216 * h ** 2 * np.sqrt(-k1 - h ** 2) * (
                                      2 * k1 ** 15 * L * np.sin(L * np.sqrt(-k1)) - k1 ** 14 * np.sin(L * np.sqrt(-k1)) * np.cosh(
                                  L * np.sqrt(-k1 - h ** 2)) + k1 ** 14 * np.sin(L * np.sqrt(-k1)) + (-k1) ** (29 / 2) * np.cos(
                                  L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2)) - (-k1) ** (29 / 2) * np.cos(
                                  L * np.sqrt(-k1)))) / (2 * (-k1) ** (3 / 2) * (-k1 - h ** 2) ** (3 / 2) * (
                        16777216 * k1 ** 13 + 62914560 * k1 ** 12 * (k1 + h ** 2) + 108003328 * k1 ** 11 * (
                        k1 + h ** 2) ** 2 + 112459776 * k1 ** 10 * (k1 + h ** 2) ** 3 + 79298560 * k1 ** 9 * (
                                k1 + h ** 2) ** 4 + 40009728 * k1 ** 8 * (
                                k1 + h ** 2) ** 5 + 14868480 * k1 ** 7 * (
                                k1 + h ** 2) ** 6 + 4122624 * k1 ** 6 * (
                                k1 + h ** 2) ** 7 + 853248 * k1 ** 5 * (k1 + h ** 2) ** 8 + 130240 * k1 ** 4 * (
                                k1 + h ** 2) ** 9 + 14256 * k1 ** 3 * (k1 + h ** 2) ** 10 + 1060 * k1 ** 2 * (
                                k1 + h ** 2) ** 11 + 48 * k1 * (k1 + h ** 2) ** 12 + (k1 + h ** 2) ** 13))
                T[3, 0, 2] = k1 * h * np.cos(L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1 - h ** 2)) / np.sqrt(-k1 - h ** 2) - h * np.sqrt(
                    -k1) * np.sin(L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1 - h ** 2))
                T[3, 0, 3] = k1 * h * np.sin(L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1 - h ** 2)) / (
                        np.sqrt(-k1) * np.sqrt(-k1 - h ** 2)) + h * np.cos(L * np.sqrt(-k1)) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - h * np.cos(L * np.sqrt(-k1))
                T[3, 1, 2] = -h * (k1 * (np.cosh(L * np.sqrt(-k1 - h ** 2)) - 1) * np.cos(L * np.sqrt(-k1)) - np.sqrt(-k1) * np.sqrt(
                    -k1 - h ** 2) * np.sin(L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1 - h ** 2))) / (k1 + h ** 2)
                T[3, 1, 3] = h * (np.sqrt(-k1) * (np.cosh(L * np.sqrt(-k1 - h ** 2)) - 1) * np.sin(L * np.sqrt(-k1)) - np.sqrt(
                    -k1 - h ** 2) * np.cos(L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1 - h ** 2))) / (k1 + h ** 2)
                T[3, 2, 4] = -(k1 * (
                        2 * L * h ** 2 * (-k1 - h ** 2) ** (3 / 2) + L * (-k1 - h ** 2) ** (5 / 2) + 2 * h ** 2 * (
                        k1 + h ** 2) * np.sinh(L * np.sqrt(-k1 - h ** 2))) * np.cos(L * np.sqrt(-k1)) + np.sqrt(-k1) * (
                                       -k1 - h ** 2) ** (3 / 2) * (
                                       k1 + 2 * h ** 2 * np.cosh(L * np.sqrt(-k1 - h ** 2)) - h ** 2) * np.sin(
                    L * np.sqrt(-k1))) / (2 * (-k1 - h ** 2) ** (5 / 2))
                T[3, 3, 4] = -(2 * k1 ** 6 * L * (k1 + h ** 2) ** 3 * np.sin(L * np.sqrt(-k1)) + 2 * k1 ** 6 * h ** 2 * (
                        k1 + h ** 2) ** 2 * np.sin(L * np.sqrt(-k1)) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - 2 * k1 ** 6 * h ** 2 * (k1 + h ** 2) ** 2 * np.sin(
                    L * np.sqrt(-k1)) + k1 ** 5 * L * (k1 + h ** 2) ** 4 * np.sin(L * np.sqrt(-k1)) / 2 + k1 ** 5 * h ** 2 * (
                                       -k1 - h ** 2) ** (5 / 2) * np.sin(L * np.sqrt(-k1)) * np.sinh(
                    L * np.sqrt(-k1 - h ** 2)) - k1 ** 5 * h ** 2 * (k1 + h ** 2) ** 3 * np.sin(
                    L * np.sqrt(-k1)) - 4 * L * h ** 2 * (-k1) ** (13 / 2) * (k1 + h ** 2) ** 2 * np.cos(
                    L * np.sqrt(-k1)) + L * h ** 2 * (-k1) ** (11 / 2) * (k1 + h ** 2) ** 3 * np.cos(
                    L * np.sqrt(-k1)) + 4 * h ** 2 * (-k1) ** (13 / 2) * (-k1 - h ** 2) ** (3 / 2) * np.cos(
                    L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1 - h ** 2)) + 2 * h ** 2 * (-k1) ** (11 / 2) * (-k1 - h ** 2) ** (
                                       5 / 2) * np.cos(L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1 - h ** 2)) - 2 * h ** 2 * (
                                   -k1) ** (11 / 2) * (k1 + h ** 2) ** 2 * np.cos(L * np.sqrt(-k1)) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) + 2 * h ** 2 * (-k1) ** (11 / 2) * (k1 + h ** 2) ** 2 * np.cos(
                    L * np.sqrt(-k1)) + h ** 2 * (-k1) ** (9 / 2) * (k1 + h ** 2) ** 3 * np.cos(L * np.sqrt(-k1)) * np.cosh(
                    L * np.sqrt(-k1 - h ** 2)) - h ** 2 * (-k1) ** (9 / 2) * (k1 + h ** 2) ** 3 * np.cos(L * np.sqrt(-k1))) / (
                                     (-k1) ** (9 / 2) * (k1 + h ** 2) ** 3 * (5 * k1 + h ** 2))
    return T