from numba import njit
import numpy as np
from numba.typed import List as nList


@njit
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
        R[0, 0] = np.cos(np.sqrt(k1) * L)
        R[0, 1] = np.sin(np.sqrt(k1) * L) / np.sqrt(k1)
        R[1, 0] = -np.sqrt(k1) * np.sin(np.sqrt(k1) * L)
        R[1, 1] = np.cos(np.sqrt(k1) * L)
        R[2, 2] = np.cosh(np.sqrt(k1) * L)
        R[2, 3] = np.sinh(np.sqrt(k1) * L) / np.sqrt(k1)
        R[3, 2] = np.sqrt(k1) * np.sinh(np.sqrt(k1) * L)
        R[3, 3] = np.cosh(np.sqrt(k1) * L)
    if k1 < 0:
        R[0, 0] = np.cosh(L * np.sqrt(-k1))
        R[0, 1] = np.sinh(L * np.sqrt(-k1)) / np.sqrt(-k1)
        R[1, 0] = np.sqrt(-k1) * np.sinh(L * np.sqrt(-k1))
        R[1, 1] = np.cosh(L * np.sqrt(-k1))
        R[2, 2] = np.cos(L * np.sqrt(-k1))
        R[2, 3] = np.sin(L * np.sqrt(-k1)) / np.sqrt(-k1)
        R[3, 2] = -np.sqrt(-k1) * np.sin(L * np.sqrt(-k1))
        R[3, 3] = np.cos(L * np.sqrt(-k1))
    return R


@njit
def compute_transport_multipole_tensor(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0]
    k1: float = element_parameters[1]
    k2: float = element_parameters[2]
    T = np.zeros((6, 6, 6))
    if k1 == 0:
        T[0, 0, 0] = -k2 * L ** 2 / 2
        T[0, 0, 1] = -k2 * L ** 3 / 3
        T[0, 1, 1] = -k2 * L ** 4 / 12
        T[0, 2, 2] = k2 * L ** 2 / 2
        T[0, 2, 3] = k2 * L ** 3 / 3
        T[0, 3, 3] = k2 * L ** 4 / 12
        T[1, 0, 0] = -k2 * L
        T[1, 0, 1] = -k2 * L ** 2
        T[1, 1, 1] = -k2 * L ** 3 / 3
        T[1, 2, 2] = k2 * L
        T[1, 2, 3] = k2 * L ** 2
        T[1, 3, 3] = k2 * L ** 3 / 3
        T[2, 0, 2] = k2 * L ** 2
        T[2, 0, 3] = k2 * L ** 3 / 3
        T[2, 1, 2] = k2 * L ** 3 / 3
        T[2, 1, 3] = k2 * L ** 4 / 6
        T[3, 0, 2] = 2 * k2 * L
        T[3, 0, 3] = k2 * L ** 2
        T[3, 1, 2] = k2 * L ** 2
        T[3, 1, 3] = 2 * k2 * L ** 3 / 3
    if k1 > 0:
        T[0, 0, 0] = k2 * (np.cos(np.sqrt(k1) * L) ** 2 + np.cos(np.sqrt(k1) * L) - 2) / (3 * k1)
        T[0, 0, 1] = 2 * k2 * (np.cos(np.sqrt(k1) * L) - 1) * np.sin(np.sqrt(k1) * L) / (3 * k1 ** (3 / 2))
        T[0, 0, 4] = np.sqrt(k1) * L * np.sin(np.sqrt(k1) * L) / 2
        T[0, 1, 1] = k2 * (np.sin(np.sqrt(k1) * L) ** 2 + 2 * np.cos(np.sqrt(k1) * L) - 2) / (3 * k1 ** 2)
        T[0, 1, 4] = -L * np.cos(np.sqrt(k1) * L) / 2 + np.sin(np.sqrt(k1) * L) / (2 * np.sqrt(k1))
        T[0, 2, 2] = k2 * (-3 * np.cos(np.sqrt(k1) * L) + np.cosh(np.sqrt(k1) * L) ** 2 + 2) / (5 * k1)
        T[0, 2, 3] = -2 * k2 * (np.sin(np.sqrt(k1) * L) - np.sinh(2 * np.sqrt(k1) * L) / 2) / (5 * k1 ** (3 / 2))
        T[0, 3, 3] = k2 * (2 * np.cos(np.sqrt(k1) * L) + np.cosh(np.sqrt(k1) * L) ** 2 - 3) / (5 * k1 ** 2)
        T[1, 0, 0] = -k2 * (np.sin(np.sqrt(k1) * L) + np.sin(2 * np.sqrt(k1) * L)) / (3 * np.sqrt(k1))
        T[1, 0, 1] = 2 * k2 * (-np.cos(np.sqrt(k1) * L) + np.cos(2 * np.sqrt(k1) * L)) / (3 * k1)
        T[1, 0, 4] = np.sqrt(k1) * np.sin(np.sqrt(k1) * L) / 2 + k1 * L * np.cos(np.sqrt(k1) * L) / 2
        T[1, 1, 1] = 2 * k2 * (np.cos(np.sqrt(k1) * L) - 1) * np.sin(np.sqrt(k1) * L) / (3 * k1 ** (3 / 2))
        T[1, 1, 4] = np.sqrt(k1) * L * np.sin(np.sqrt(k1) * L) / 2
        T[1, 2, 2] = k2 * (3 * np.sin(np.sqrt(k1) * L) + np.sinh(2 * np.sqrt(k1) * L)) / (5 * np.sqrt(k1))
        T[1, 2, 3] = -2 * k2 * (np.cos(np.sqrt(k1) * L) - np.cosh(2 * np.sqrt(k1) * L)) / (5 * k1)
        T[1, 3, 3] = -2 * k2 * (np.sin(np.sqrt(k1) * L) - np.sinh(2 * np.sqrt(k1) * L) / 2) / (5 * k1 ** (3 / 2))
        T[2, 0, 2] = 2 * k2 * (
                2 * np.sin(np.sqrt(k1) * L) * np.sinh(np.sqrt(k1) * L) - np.cos(np.sqrt(k1) * L) * np.cosh(
            np.sqrt(k1) * L) + np.cosh(np.sqrt(k1) * L)) / (5 * k1)
        T[2, 0, 3] = -2 * k2 * (
                -2 * np.sin(np.sqrt(k1) * L) * np.cosh(np.sqrt(k1) * L) + np.cos(np.sqrt(k1) * L) * np.sinh(
            np.sqrt(k1) * L) + np.sinh(np.sqrt(k1) * L)) / (5 * k1 ** (3 / 2))
        T[2, 1, 2] = -2 * k2 * (
                np.sin(np.sqrt(k1) * L) * np.cosh(np.sqrt(k1) * L) + 2 * np.cos(np.sqrt(k1) * L) * np.sinh(
            np.sqrt(k1) * L) - 3 * np.sinh(np.sqrt(k1) * L)) / (5 * k1 ** (3 / 2))
        T[2, 1, 3] = -2 * k2 * (
                np.sin(np.sqrt(k1) * L) * np.sinh(np.sqrt(k1) * L) + 2 * np.cos(np.sqrt(k1) * L) * np.cosh(
            np.sqrt(k1) * L) - 2 * np.cosh(np.sqrt(k1) * L)) / (5 * k1 ** 2)
        T[2, 2, 4] = -np.sqrt(k1) * L * np.sinh(np.sqrt(k1) * L) / 2
        T[2, 3, 4] = -L * np.cosh(np.sqrt(k1) * L) / 2 + np.sinh(np.sqrt(k1) * L) / (2 * np.sqrt(k1))
        T[3, 0, 2] = 2 * k2 * (
                3 * np.sin(np.sqrt(k1) * L) * np.cosh(np.sqrt(k1) * L) + np.cos(np.sqrt(k1) * L) * np.sinh(
            np.sqrt(k1) * L) + np.sinh(np.sqrt(k1) * L)) / (5 * np.sqrt(k1))
        T[3, 0, 3] = 2 * k2 * (
                3 * np.sin(np.sqrt(k1) * L) * np.sinh(np.sqrt(k1) * L) + np.cos(np.sqrt(k1) * L) * np.cosh(
            np.sqrt(k1) * L) - np.cosh(np.sqrt(k1) * L)) / (5 * k1)
        T[3, 1, 2] = 2 * k2 * (
                np.sin(np.sqrt(k1) * L) * np.sinh(np.sqrt(k1) * L) - 3 * np.cos(np.sqrt(k1) * L) * np.cosh(
            np.sqrt(k1) * L) + 3 * np.cosh(np.sqrt(k1) * L)) / (5 * k1)
        T[3, 1, 3] = 2 * k2 * (
                np.sin(np.sqrt(k1) * L) * np.cosh(np.sqrt(k1) * L) - 3 * np.cos(np.sqrt(k1) * L) * np.sinh(
            np.sqrt(k1) * L) + 2 * np.sinh(np.sqrt(k1) * L)) / (5 * k1 ** (3 / 2))
        T[3, 2, 4] = -np.sqrt(k1) * np.sinh(np.sqrt(k1) * L) / 2 - k1 * L * np.cosh(np.sqrt(k1) * L) / 2
        T[3, 3, 4] = -np.sqrt(k1) * L * np.sinh(np.sqrt(k1) * L) / 2
    if k1 < 0:
        T[0, 0, 0] = k2 * (np.cos(np.sqrt(k1) * L) ** 2 + np.cosh(L * np.sqrt(-k1)) - 2) / (3 * k1)
        T[0, 0, 1] = -2 * k2 * (1 - np.cosh(L * np.sqrt(-k1))) * np.sinh(L * np.sqrt(-k1)) / (3 * (-k1) ** (3 / 2))
        T[0, 0, 4] = -L * np.sqrt(-k1) * np.sinh(L * np.sqrt(-k1)) / 2
        T[0, 1, 1] = k2 * (-np.sinh(L * np.sqrt(-k1)) ** 2 + 2 * np.cosh(L * np.sqrt(-k1)) - 2) / (3 * k1 ** 2)
        T[0, 1, 4] = -(-L * np.sqrt(-k1) * np.cosh(L * np.sqrt(-k1)) + np.sinh(L * np.sqrt(-k1))) / (2 * np.sqrt(-k1))
        T[0, 2, 2] = k2 * (np.cosh(np.sqrt(k1) * L) ** 2 - 3 * np.cosh(L * np.sqrt(-k1)) + 2) / (5 * k1)
        T[0, 2, 3] = -2 * k2 * (-np.sin(2 * L * np.sqrt(-k1)) / 2 + np.sinh(L * np.sqrt(-k1))) / (5 * (-k1) ** (3 / 2))
        T[0, 3, 3] = k2 * (np.cosh(np.sqrt(k1) * L) ** 2 + 2 * np.cosh(L * np.sqrt(-k1)) - 3) / (5 * k1 ** 2)
        T[1, 0, 0] = -k2 * (
                    -np.sqrt(-k1) * np.sinh(L * np.sqrt(-k1)) - np.sqrt(-k1) * np.sinh(2 * L * np.sqrt(-k1))) / (3 * k1)
        T[1, 0, 1] = 2 * k2 * (np.cosh(L * np.sqrt(-k1)) - np.cosh(2 * L * np.sqrt(-k1))) / (3 * k1)
        T[1, 0, 4] = k1 * L * np.cosh(L * np.sqrt(-k1)) / 2 - np.sqrt(-k1) * np.sinh(L * np.sqrt(-k1)) / 2
        T[1, 1, 1] = -2 * k2 * (np.cosh(L * np.sqrt(-k1)) - 1) * np.sinh(L * np.sqrt(-k1)) / (3 * (-k1) ** (3 / 2))
        T[1, 1, 4] = L * np.sqrt(-k1) * np.sinh(L * np.sqrt(-k1)) / 2
        T[1, 2, 2] = k2 * (
                    -np.sqrt(-k1) * np.sin(2 * L * np.sqrt(-k1)) - 3 * np.sqrt(-k1) * np.sinh(L * np.sqrt(-k1))) / (
                                 5 * k1)
        T[1, 2, 3] = -2 * k2 * (np.cos(2 * L * np.sqrt(-k1)) - np.cosh(L * np.sqrt(-k1))) / (5 * k1)
        T[1, 3, 3] = 2 * k2 * (
                    -np.sqrt(-k1) * np.sin(2 * L * np.sqrt(-k1)) / 2 + np.sqrt(-k1) * np.sinh(L * np.sqrt(-k1))) / (
                             5 * k1 ** 2)
        T[2, 0, 2] = 2 * k2 * (
                -2 * np.sin(L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1)) - np.cos(L * np.sqrt(-k1)) * np.cosh(
            L * np.sqrt(-k1)) + np.cos(
            L * np.sqrt(-k1))) / (5 * k1)
        T[2, 0, 3] = -2 * k2 * (
                np.sin(L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1)) + np.sin(L * np.sqrt(-k1)) - 2 * np.cos(
            L * np.sqrt(-k1)) * np.sinh(
            L * np.sqrt(-k1))) / (5 * (-k1) ** (3 / 2))
        T[2, 1, 2] = -2 * k2 * (
                2 * np.sin(L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1)) - 3 * np.sin(L * np.sqrt(-k1)) + np.cos(
            L * np.sqrt(-k1)) * np.sinh(
            L * np.sqrt(-k1))) / (5 * (-k1) ** (3 / 2))
        T[2, 1, 3] = -2 * k2 * (
                -np.sin(L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1)) + 2 * np.cos(L * np.sqrt(-k1)) * np.cosh(
            L * np.sqrt(-k1)) - 2 * np.cos(
            L * np.sqrt(-k1))) / (5 * k1 ** 2)
        T[2, 2, 4] = L * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) / 2
        T[2, 3, 4] = L * np.cos(L * np.sqrt(-k1)) / 2 - np.sin(L * np.sqrt(-k1)) / (2 * np.sqrt(-k1))
        T[3, 0, 2] = 2 * k2 * (
                np.sin(L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1)) + np.sin(L * np.sqrt(-k1)) + 3 * np.cos(
            L * np.sqrt(-k1)) * np.sinh(
            L * np.sqrt(-k1))) / (5 * np.sqrt(-k1))
        T[3, 0, 3] = 2 * k2 * (
                3 * np.sin(L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1)) - np.cos(L * np.sqrt(-k1)) * np.cosh(
            L * np.sqrt(-k1)) + np.cos(
            L * np.sqrt(-k1))) / (5 * k1)
        T[3, 1, 2] = 2 * k2 * (
                np.sin(L * np.sqrt(-k1)) * np.sinh(L * np.sqrt(-k1)) + 3 * np.cos(L * np.sqrt(-k1)) * np.cosh(
            L * np.sqrt(-k1)) - 3 * np.cos(L * np.sqrt(-k1))) / (5 * k1)
        T[3, 1, 3] = -2 * k2 * (
                -3 * np.sin(L * np.sqrt(-k1)) * np.cosh(L * np.sqrt(-k1)) + 2 * np.sin(L * np.sqrt(-k1)) + np.cos(
            L * np.sqrt(-k1)) * np.sinh(
            L * np.sqrt(-k1))) / (5 * (-k1) ** (3 / 2))
        T[3, 2, 4] = -k1 * L * np.cos(L * np.sqrt(-k1)) / 2 + np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) / 2
        T[3, 3, 4] = -L * np.sqrt(-k1) * np.sin(L * np.sqrt(-k1)) / 2
    return T
