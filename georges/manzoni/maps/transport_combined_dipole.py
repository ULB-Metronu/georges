from numba import njit
import numpy as np
from numpy import cos, sin, cosh, sinh, sqrt


@njit(cache=True)
def compute_transport_combined_dipole_matrix(L: float, alpha: float, K1: float, K2: float, *args) -> np.ndarray:
    h = alpha / L
    kx2 = h ** 2 + K1
    ky2 = -K1
    R = np.zeros((6, 6))

    if h ** 2 + K1 > 0.0 and K1 > 0 and h > 1e-6:
        kx = sqrt(kx2)
        ky = sqrt(-ky2)

        R[0, 0] = cos(L * (kx))
        R[0, 1] = (kx ** (-1)) * sin(L * (kx))
        R[0, 5] = h * (1 - cos(L * (kx))) * (kx ** (-2))
        R[1, 0] = (kx ** (-1)) * (-kx ** 2) * sin(L * (kx))
        R[1, 1] = cos(L * (kx))
        R[1, 5] = h * (kx ** (-1)) * sin(L * (kx))
        R[2, 2] = cosh(sqrt(K1) * L)
        R[2, 3] = sinh(sqrt(K1) * L) / sqrt(K1)
        R[3, 2] = sqrt(K1) * sinh(sqrt(K1) * L)
        R[3, 3] = cosh(sqrt(K1) * L)
        R[4, 4] = 1
        R[5, 5] = 1
        return R

    if h ** 2 + K1 > 0.0 and K1 < 0 and h > 1e-6:
        kx = sqrt(kx2)
        ky = sqrt(ky2)

        R[0, 0] = cos(L * (kx))
        R[0, 1] = (kx ** (-1)) * sin(L * (kx))
        R[0, 5] = h * (1 - cos(L * (kx))) * (kx ** (-2))
        R[1, 0] = (kx ** (-1)) * (-kx ** 2) * sin(L * (kx))
        R[1, 1] = cos(L * (kx))
        R[1, 5] = h * (kx ** (-1)) * sin(L * (kx))
        R[2, 2] = cos(sqrt(-K1) * L)
        R[2, 3] = sin(sqrt(-K1) * L) / sqrt(-K1)
        R[3, 2] = (K1 * sin(sqrt(-K1) * L)) / sqrt(-K1)
        R[3, 3] = cos(sqrt(-K1) * L)
        R[4, 4] = 1
        R[5, 5] = 1
        return R

    if h ** 2 + K1 < 0.0 and h > 1e-6:
        kx = sqrt(-kx2)
        ky = sqrt(ky2)

        R[0, 0] = cosh(L * (kx))
        R[0, 1] = (kx ** (-1)) * sinh(L * (kx))
        R[0, 5] = h * (1 - cosh(L * (kx))) * (kx ** (-2))
        R[1, 0] = (kx) * sinh(L * (kx))
        R[1, 1] = cosh(L * (kx))
        R[1, 5] = h * (kx ** (-1)) * sinh(L * (kx))
        R[2, 2] = cos(sqrt(-K1) * L)
        R[2, 3] = sin(sqrt(-K1) * L) / sqrt(-K1)
        R[3, 2] = (K1 * sin(sqrt(-K1) * L)) / sqrt(-K1)
        R[3, 3] = cos(sqrt(-K1) * L)
        R[4, 4] = 1
        R[5, 5] = 1
        return R

    if np.abs(kx2 - 4 * ky2) < 1e-3:
        R[0, 0] = 1
        R[0, 1] = L
        R[0, 5] = (h * L ** 2) / 2
        R[1, 0] = 0
        R[1, 1] = 1
        R[1, 5] = h * L
        R[2, 2] = 1 - (h ** 2 * L ** 2) / 2 + (h ** 4 * L ** 4) / 24 - (h ** 6 * L ** 6) / 720
        R[2, 3] = L - (h ** 2 * L ** 3) / 6 + (h ** 4 * L ** 5) / 120 - (h ** 6 * L ** 7) / 5040
        R[3, 2] = -(h ** 2 * L) + (h ** 4 * L ** 3) / 6 - (h ** 6 * L ** 5) / 120
        R[3, 3] = 1 - (h ** 2 * L ** 2) / 2 + (h ** 4 * L ** 4) / 24 - (h ** 6 * L ** 6) / 720
        R[4, 4] = 1
        R[5, 5] = 1
        return R

    if np.abs(kx2 * ky2) < 1e-3:
        R[0, 0] = 1 - (h ** 2 * L ** 2) / 2 + (h ** 4 * L ** 4) / 24 - (h ** 6 * L ** 6) / 720
        R[0, 1] = L - (h ** 2 * L ** 3) / 6 + (h ** 4 * L ** 5) / 120 - (h ** 6 * L ** 7) / 5040
        R[0, 5] = (h * L ** 2) / 2 - (h ** 3 * L ** 4) / 24 + (h ** 5 * L ** 6) / 720
        R[1, 0] = -(h ** 2 * L) + (h ** 4 * L ** 3) / 6 - (h ** 6 * L ** 5) / 120
        R[1, 1] = 1 - (h ** 2 * L ** 2) / 2 + (h ** 4 * L ** 4) / 24 - (h ** 6 * L ** 6) / 720
        R[1, 5] = h * L - (h ** 3 * L ** 3) / 6 + (h ** 5 * L ** 5) / 120
        R[2, 2] = 1
        R[2, 3] = L
        R[3, 2] = 0
        R[3, 3] = 1
        R[4, 4] = 1
        R[5, 5] = 1
        return R

    if h < 1e-6:

        R[0, 0] = 1 - (K1 * L ** 2) / 2 + (K1 ** 2 * L ** 4) / 24 - (K1 ** 3 * L ** 6) / 720 + (
                K1 ** 4 * L ** 8) / 40320 - (K1 ** 5 * L ** 10) / 3628800 + (K1 ** 6 * L ** 12) / 479001600
        R[0, 1] = L - (K1 * L ** 3) / 6 + (K1 ** 2 * L ** 5) / 120 - (K1 ** 3 * L ** 7) / 5040 + (
                K1 ** 4 * L ** 9) / 362880 - (K1 ** 5 * L ** 11) / 39916800 + (K1 ** 6 * L ** 13) / 6227020800
        R[0, 5] = 0
        R[1, 0] = -(K1 * L) + (K1 ** 2 * L ** 3) / 6 - (K1 ** 3 * L ** 5) / 120 + (K1 ** 4 * L ** 7) / 5040 - (
                K1 ** 5 * L ** 9) / 362880 + (K1 ** 6 * L ** 11) / 39916800
        R[1, 1] = 1 - (K1 * L ** 2) / 2 + (K1 ** 2 * L ** 4) / 24 - (K1 ** 3 * L ** 6) / 720 + (
                K1 ** 4 * L ** 8) / 40320 - (K1 ** 5 * L ** 10) / 3628800 + (K1 ** 6 * L ** 12) / 479001600
        R[1, 5] = 0
        R[2, 2] = 1 + (K1 * L ** 2) / 2 + (K1 ** 2 * L ** 4) / 24 + (K1 ** 3 * L ** 6) / 720 + (
                K1 ** 4 * L ** 8) / 40320 + (K1 ** 5 * L ** 10) / 3628800 + (K1 ** 6 * L ** 12) / 479001600
        R[2, 3] = L + (K1 * L ** 3) / 6 + (K1 ** 2 * L ** 5) / 120 + (K1 ** 3 * L ** 7) / 5040 + (
                K1 ** 4 * L ** 9) / 362880 + (K1 ** 5 * L ** 11) / 39916800 + (K1 ** 6 * L ** 13) / 6227020800
        R[3, 2] = K1 * L + (K1 ** 2 * L ** 3) / 6 + (K1 ** 3 * L ** 5) / 120 + (K1 ** 4 * L ** 7) / 5040 + (
                K1 ** 5 * L ** 9) / 362880 + (K1 ** 6 * L ** 11) / 39916800
        R[3, 3] = 1 + (K1 * L ** 2) / 2 + (K1 ** 2 * L ** 4) / 24 + (K1 ** 3 * L ** 6) / 720 + (
                K1 ** 4 * L ** 8) / 40320 + (K1 ** 5 * L ** 10) / 3628800 + (K1 ** 6 * L ** 12) / 479001600
        R[4, 4] = 1
        R[5, 5] = 1
        return R


@njit(cache=True)
def compute_transport_combined_dipole_tensor(L: float, alpha: float, K1: float, K2: float, *args) -> np.ndarray:
    h = alpha / L
    kx2 = h ** 2 + K1
    ky2 = -K1
    T = np.zeros((6, 6, 6))

    if h ** 2 + K1 > 0.0 and K1 > 0 and h > 1e-6:
        kx = sqrt(kx2)
        ky = sqrt(-ky2)

        T[0, 0, 0] = (h ** 3 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                    (1 - cos(L * (kx))) * (kx ** (-2)) + (kx ** (-2)) * sin(L * (kx)) ** 2)) / 3 + (h * (kx ** 2) * (
                    (1 - cos(L * (kx))) * (kx ** (-2)) + (
                        -((1 - cos(L * (kx))) * (kx ** (-2))) - (kx ** (-2)) * sin(L * (kx)) ** 2) / 3)) / 2
        T[0, 0, 1] = h * (kx ** (-1)) * sin(L * (kx)) - (h * (1 - cos(L * (kx))) * (kx ** (-1)) * sin(L * (kx))) / 3 + (
                    2 * h ** 3 * (1 - cos(L * (kx))) * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * sin(
                L * (kx))) / (3 * (kx ** 3))
        T[0, 0, 5] = (h ** 2 * L * (kx ** (-1)) * (1 + (kx ** 2) / h ** 2) * sin(L * (kx))) / 2 - h ** 2 * (
                    (1 - cos(L * (kx))) * (kx ** (-2)) + (
                        -((1 - cos(L * (kx))) * (kx ** (-2))) - (kx ** (-2)) * sin(L * (kx)) ** 2) / 3) + 2 * h ** 4 * (
                                 kx ** (-2)) * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 (L * (kx ** (-1)) * sin(L * (kx))) / 2 + (
                                     -((1 - cos(L * (kx))) * (kx ** (-2))) - (kx ** (-2)) * sin(L * (kx)) ** 2) / 3)
        T[0, 1, 1] = (h * ((1 - cos(L * (kx))) * (kx ** (-2)) + (kx ** (-2)) * sin(L * (kx)) ** 2)) / 6 + h ** 3 * (
                    kx ** (-2)) * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 (1 - cos(L * (kx))) * (kx ** (-2)) + (
                                     -((1 - cos(L * (kx))) * (kx ** (-2))) - (kx ** (-2)) * sin(L * (kx)) ** 2) / 3)
        T[0, 1, 5] = (h ** 2 * (1 - cos(L * (kx))) * sin(L * (kx))) / (3 * (kx ** 3)) + (
                    h ** 2 * (kx ** (-2)) * (1 + (kx ** 2) / h ** 2) * (
                        -(L * cos(L * (kx))) + (kx ** (-1)) * sin(L * (kx)))) / 2 + 2 * h ** 4 * (kx ** (-2)) * (
                                 -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 -((1 - cos(L * (kx))) * sin(L * (kx))) / (3 * (kx ** 3)) + (
                                     (kx ** (-2)) * (-(L * cos(L * (kx))) + (kx ** (-1)) * sin(L * (kx)))) / 2)
        T[0, 5, 5] = -(h * (1 - cos(L * (kx))) * (kx ** (-2))) + h ** 3 * (kx ** (-2)) * (1 + (kx ** 2) / h ** 2) * (
                    (1 - cos(L * (kx))) * (kx ** (-2)) - (L * (kx ** (-1)) * sin(L * (kx))) / 2) + (
                                 h ** 5 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                     (4 * (1 - cos(L * (kx))) * (kx ** (-2))) / 3 - L * (kx ** (-1)) * sin(L * (kx)) + (
                                         (kx ** (-2)) * sin(L * (kx)) ** 2) / 3)) / (kx ** 2) ** 2 + (
                                 h ** 3 * (kx ** (-2)) * ((1 - cos(L * (kx))) * (kx ** (-2)) + (
                                     -((1 - cos(L * (kx))) * (kx ** (-2))) - (kx ** (-2)) * sin(
                                 L * (kx)) ** 2) / 3)) / 2
        T[0, 2, 2] = (h * K1 * (1 - cos(L * (kx))) * (kx ** (-2))) / 2 + K2 * ((1 - cos(L * (kx))) * (kx ** (-2)) + (
                    K1 * (-2 * (1 - cos(L * (kx))) * (kx ** (-2)) + sinh(sqrt(K1) * L) ** 2 / K1)) / (h ** 2 + 5 * K1))
        T[0, 2, 3] = (2 * K2 * (
                    -((kx ** (-1)) * sin(L * (kx))) + (cosh(sqrt(K1) * L) * sinh(sqrt(K1) * L)) / sqrt(K1))) / (
                                 h ** 2 + 5 * K1)
        T[0, 3, 3] = -(h * (1 - cos(L * (kx))) * (kx ** (-2))) / 2 + (
                    K2 * (-2 * (1 - cos(L * (kx))) * (kx ** (-2)) + sinh(sqrt(K1) * L) ** 2 / K1)) / (h ** 2 + 5 * K1)
        T[1, 0, 0] = (h * (1 - cos(L * (kx))) * (kx) * sin(L * (kx))) / 3 + h * cos(L * (kx)) * (kx) * sin(L * (kx)) + (
                    h ** 3 * (1 + 2 * cos(L * (kx))) * (kx ** (-1)) * (
                        -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * sin(L * (kx))) / 3
        T[1, 0, 1] = h * cos(L * (kx)) - h * (cos(L * (kx)) ** 2 - sin(L * (kx)) ** 2) - (h * (kx ** 2) * (
                    -((1 - cos(L * (kx))) * (kx ** (-2))) + 2 * (kx ** (-2)) * sin(L * (kx)) ** 2)) / 3 + (
                                 2 * h ** 3 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                     -((1 - cos(L * (kx))) * (kx ** (-2))) + 2 * (kx ** (-2)) * sin(L * (kx)) ** 2)) / 3
        T[1, 0, 5] = (-2 * h ** 2 * (1 - cos(L * (kx))) * (kx ** (-1)) * sin(L * (kx))) / 3 + (
                    h ** 2 * (1 + (kx ** 2) / h ** 2) * (
                        L * cos(L * (kx)) + (kx ** (-1)) * sin(L * (kx)))) / 2 + 2 * h ** 4 * (kx ** (-2)) * (
                                 -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 (L * cos(L * (kx))) / 2 + ((kx ** (-1)) * sin(L * (kx))) / 6 - (
                                     2 * cos(L * (kx)) * (kx ** (-1)) * sin(L * (kx))) / 3) - h * (
                                 -(h * (1 - cos(L * (kx))) * (kx ** (-1)) * sin(L * (kx))) + h * cos(L * (kx)) * (
                                     kx ** (-1)) * sin(L * (kx)))
        T[1, 1, 1] = -(h * cos(L * (kx)) * (kx ** (-1)) * sin(L * (kx))) + (
                    h * (1 + 2 * cos(L * (kx))) * (kx ** (-1)) * sin(L * (kx))) / 6 + (
                                 2 * h ** 3 * (1 - cos(L * (kx))) * (
                                     -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * sin(L * (kx))) / (3 * (kx ** 3))
        T[1, 1, 5] = (h ** 2 * L * (kx ** (-1)) * (1 + (kx ** 2) / h ** 2) * sin(L * (kx))) / 2 + 2 * h ** 4 * (
                    kx ** (-2)) * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 ((1 - cos(L * (kx))) * (kx ** (-2))) / 3 + (L * (kx ** (-1)) * sin(L * (kx))) / 2 - (
                                     2 * (kx ** (-2)) * sin(L * (kx)) ** 2) / 3) + (h ** 2 * (
                    -((1 - cos(L * (kx))) * (kx ** (-2))) + 2 * (kx ** (-2)) * sin(L * (kx)) ** 2)) / 3 - h * (
                                 h * (1 - cos(L * (kx))) * cos(L * (kx)) * (kx ** (-2)) + h * (kx ** (-2)) * sin(
                             L * (kx)) ** 2)
        T[1, 5, 5] = -(h * (kx ** (-1)) * sin(L * (kx))) - (2 * h ** 3 * (1 - cos(L * (kx))) * sin(L * (kx))) / (
                    3 * (kx ** 3)) + (h ** 3 * (kx ** (-2)) * (1 + (kx ** 2) / h ** 2) * (
                    -(L * cos(L * (kx))) + (kx ** (-1)) * sin(L * (kx)))) / 2 + (
                                 h ** 5 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                     -(L * cos(L * (kx))) + ((kx ** (-1)) * sin(L * (kx))) / 3 + (
                                         2 * cos(L * (kx)) * (kx ** (-1)) * sin(L * (kx))) / 3)) / (kx ** 2) ** 2
        T[1, 2, 2] = (h * K1 * (kx ** (-1)) * sin(L * (kx))) / 2 + K2 * ((kx ** (-1)) * sin(L * (kx)) + (2 * K1 * (
                    -((kx ** (-1)) * sin(L * (kx))) + (cosh(sqrt(K1) * L) * sinh(sqrt(K1) * L)) / sqrt(K1))) / (
                                                                                     h ** 2 + 5 * K1))
        T[1, 2, 3] = (2 * K2 * (1 - cos(L * (kx)) + 2 * sinh(sqrt(K1) * L) ** 2)) / (h ** 2 + 5 * K1)
        T[1, 3, 3] = -(h * (kx ** (-1)) * sin(L * (kx))) / 2 + (2 * K2 * (
                    -((kx ** (-1)) * sin(L * (kx))) + (cosh(sqrt(K1) * L) * sinh(sqrt(K1) * L)) / sqrt(K1))) / (
                                 h ** 2 + 5 * K1)
        T[2, 0, 2] = (2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                    (1 - cos(L * (kx))) * cosh(sqrt(K1) * L) + 2 * sqrt(K1 * (kx ** (-2))) * sin(L * (kx)) * sinh(
                sqrt(K1) * L))) / (h ** 2 + 5 * K1) - (h * K1 * (kx ** 2) * (
                    2 * (1 - cos(L * (kx))) * cosh(sqrt(K1) * L) * (kx ** (-2)) - ((sqrt(K1) * kx) ** (-1)) * sin(
                L * (kx)) * sinh(sqrt(K1) * L))) / (h ** 2 + 5 * K1)
        T[2, 0, 3] = (h * sinh(sqrt(K1) * L)) / sqrt(K1) + (2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                    2 * cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx)) - (
                        (1 + cos(L * (kx))) * sinh(sqrt(K1) * L)) / sqrt(K1))) / (h ** 2 + 5 * K1) - h * (kx ** 2) * (
                                 ((kx ** (-2)) * sinh(sqrt(K1) * L)) / sqrt(K1) + (
                                     -(cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx))) - 2 * sqrt(K1) * (
                                         1 + cos(L * (kx))) * (kx ** (-2)) * sinh(sqrt(K1) * L)) / (h ** 2 + 5 * K1))
        T[2, 1, 2] = (h * K1 * (2 * cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx)) - (
                    (1 + cos(L * (kx))) * sinh(sqrt(K1) * L)) / sqrt(K1))) / (h ** 2 + 5 * K1) + 2 * h ** 3 * (
                                 K1 / h ** 2 + K2 / h ** 3) * (((kx ** (-2)) * sinh(sqrt(K1) * L)) / sqrt(K1) + (
                    -(cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx))) - 2 * sqrt(K1) * (1 + cos(L * (kx))) * (
                        kx ** (-2)) * sinh(sqrt(K1) * L)) / (h ** 2 + 5 * K1))
        T[2, 1, 3] = (h * (
                    (1 - cos(L * (kx))) * cosh(sqrt(K1) * L) + 2 * sqrt(K1 * (kx ** (-2))) * sin(L * (kx)) * sinh(
                sqrt(K1) * L))) / (h ** 2 + 5 * K1) + (2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                    2 * (1 - cos(L * (kx))) * cosh(sqrt(K1) * L) * (kx ** (-2)) - ((sqrt(K1) * kx) ** (-1)) * sin(
                L * (kx)) * sinh(sqrt(K1) * L))) / (h ** 2 + 5 * K1)
        T[2, 2, 5] = -(sqrt(K1) * L * sinh(sqrt(K1) * L)) / 2 + (h ** 2 * K1 * (
                    2 * (1 - cos(L * (kx))) * cosh(sqrt(K1) * L) * (kx ** (-2)) - ((sqrt(K1) * kx) ** (-1)) * sin(
                L * (kx)) * sinh(sqrt(K1) * L))) / (h ** 2 + 5 * K1) + 2 * h ** 4 * (kx ** (-2)) * (
                                 -1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * ((L * sinh(sqrt(K1) * L)) / (2 * sqrt(K1)) - (
                    (1 - cos(L * (kx))) * cosh(sqrt(K1) * L) + 2 * sqrt(K1 * (kx ** (-2))) * sin(L * (kx)) * sinh(
                sqrt(K1) * L)) / (h ** 2 + 5 * K1))
        T[2, 3, 5] = (-(L * cosh(sqrt(K1) * L)) + sinh(sqrt(K1) * L) / sqrt(K1)) / 2 + 2 * h ** 4 * (kx ** (-2)) * (
                    -1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                                 -(-(L * cosh(sqrt(K1) * L)) + sinh(sqrt(K1) * L) / sqrt(K1)) / (2 * K1) - (
                                     2 * cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx)) - (
                                         (1 + cos(L * (kx))) * sinh(sqrt(K1) * L)) / sqrt(K1)) / (
                                             h ** 2 + 5 * K1)) + h ** 2 * (
                                 ((kx ** (-2)) * sinh(sqrt(K1) * L)) / sqrt(K1) + (
                                     -(cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx))) - 2 * sqrt(K1) * (
                                         1 + cos(L * (kx))) * (kx ** (-2)) * sinh(sqrt(K1) * L)) / (h ** 2 + 5 * K1))
        T[3, 0, 2] = -(h * sqrt(K1) * cos(L * (kx)) * sinh(sqrt(K1) * L)) + (
                    2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                        (h ** 2 + 3 * K1) * cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx)) + sqrt(K1) * (
                            1 + cos(L * (kx))) * sinh(sqrt(K1) * L))) / (h ** 2 + 5 * K1) - (h * K1 * (kx ** 2) * (
                    cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx)) - (cos(L * (kx)) * sinh(sqrt(K1) * L)) / sqrt(
                K1) + 2 * sqrt(K1) * (1 - cos(L * (kx))) * (kx ** (-2)) * sinh(sqrt(K1) * L))) / (h ** 2 + 5 * K1)
        T[3, 0, 3] = h * cosh(sqrt(K1) * L) - h * cos(L * (kx)) * cosh(sqrt(K1) * L) + (
                    2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                        -((1 - cos(L * (kx))) * cosh(sqrt(K1) * L)) + (h ** 2 + 3 * K1) * (
                            (sqrt(K1) * kx) ** (-1)) * sin(L * (kx)) * sinh(sqrt(K1) * L))) / (h ** 2 + 5 * K1) - h * (
                                 kx ** 2) * (cosh(sqrt(K1) * L) * (kx ** (-2)) + (
                    -(cos(L * (kx)) * cosh(sqrt(K1) * L)) - 2 * K1 * (1 + cos(L * (kx))) * cosh(sqrt(K1) * L) * (
                        kx ** (-2)) + sqrt(K1 * (kx ** (-2))) * sin(L * (kx)) * sinh(sqrt(K1) * L)) / (h ** 2 + 5 * K1))
        T[3, 1, 2] = -(h * sqrt(K1 * (kx ** (-2))) * sin(L * (kx)) * sinh(sqrt(K1) * L)) + (h * K1 * (
                    -((1 - cos(L * (kx))) * cosh(sqrt(K1) * L)) + (h ** 2 + 3 * K1) * ((sqrt(K1) * kx) ** (-1)) * sin(
                L * (kx)) * sinh(sqrt(K1) * L))) / (h ** 2 + 5 * K1) + 2 * h ** 3 * (K1 / h ** 2 + K2 / h ** 3) * (
                                 cosh(sqrt(K1) * L) * (kx ** (-2)) + (
                                     -(cos(L * (kx)) * cosh(sqrt(K1) * L)) - 2 * K1 * (1 + cos(L * (kx))) * cosh(
                                 sqrt(K1) * L) * (kx ** (-2)) + sqrt(K1 * (kx ** (-2))) * sin(L * (kx)) * sinh(
                                 sqrt(K1) * L)) / (h ** 2 + 5 * K1))
        T[3, 1, 3] = -(h * cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx))) + (h * (
                    (h ** 2 + 3 * K1) * cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx)) + sqrt(K1) * (
                        1 + cos(L * (kx))) * sinh(sqrt(K1) * L))) / (h ** 2 + 5 * K1) + (
                                 2 * h ** 3 * (K1 / h ** 2 + K2 / h ** 3) * (
                                     cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx)) - (
                                         cos(L * (kx)) * sinh(sqrt(K1) * L)) / sqrt(K1) + 2 * sqrt(K1) * (
                                                 1 - cos(L * (kx))) * (kx ** (-2)) * sinh(sqrt(K1) * L))) / (
                                 h ** 2 + 5 * K1)
        T[3, 2, 5] = -(h ** 2 * sqrt(K1) * (1 - cos(L * (kx))) * (kx ** (-2)) * sinh(sqrt(K1) * L)) - (
                    K1 * (L * cosh(sqrt(K1) * L) + sinh(sqrt(K1) * L) / sqrt(K1))) / 2 + (h ** 2 * K1 * (
                    cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx)) - (cos(L * (kx)) * sinh(sqrt(K1) * L)) / sqrt(
                K1) + 2 * sqrt(K1) * (1 - cos(L * (kx))) * (kx ** (-2)) * sinh(sqrt(K1) * L))) / (
                                 h ** 2 + 5 * K1) + 2 * h ** 4 * (K1 / h ** 2 + K2 / h ** 3) * (kx ** (-2)) * (
                                 (L * cosh(sqrt(K1) * L)) / 2 + sinh(sqrt(K1) * L) / (2 * sqrt(K1)) + (
                                     -((h ** 2 + 3 * K1) * cosh(sqrt(K1) * L) * (kx ** (-1)) * sin(L * (kx))) - sqrt(
                                 K1) * (1 + cos(L * (kx))) * sinh(sqrt(K1) * L)) / (h ** 2 + 5 * K1))
        T[3, 3, 5] = -(h ** 2 * (1 - cos(L * (kx))) * cosh(sqrt(K1) * L) * (kx ** (-2))) - (
                    sqrt(K1) * L * sinh(sqrt(K1) * L)) / 2 + h ** 2 * (cosh(sqrt(K1) * L) * (kx ** (-2)) + (
                    -(cos(L * (kx)) * cosh(sqrt(K1) * L)) - 2 * K1 * (1 + cos(L * (kx))) * cosh(sqrt(K1) * L) * (
                        kx ** (-2)) + sqrt(K1 * (kx ** (-2))) * sin(L * (kx)) * sinh(sqrt(K1) * L)) / (
                                                                                   h ** 2 + 5 * K1)) + 2 * h ** 4 * (
                                 K1 / h ** 2 + K2 / h ** 3) * (kx ** (-2)) * (
                                 (L * sinh(sqrt(K1) * L)) / (2 * sqrt(K1)) - (
                                     -((1 - cos(L * (kx))) * cosh(sqrt(K1) * L)) + (h ** 2 + 3 * K1) * (
                                         (sqrt(K1) * kx) ** (-1)) * sin(L * (kx)) * sinh(sqrt(K1) * L)) / (
                                             h ** 2 + 5 * K1))
        return T

    if h ** 2 + K1 > 0.0 and K1 < 0 and h > 1e-6:
        kx = sqrt(kx2)
        ky = sqrt(ky2)

        T[0, 0, 0] = (h ** 3 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                    (1 - cos(L * (kx))) * (kx ** (-2)) + (kx ** (-2)) * sin(L * (kx)) ** 2)) / 3 + (h * (kx ** 2) * (
                    (1 - cos(L * (kx))) * (kx ** (-2)) + (
                        -((1 - cos(L * (kx))) * (kx ** (-2))) - (kx ** (-2)) * sin(L * (kx)) ** 2) / 3)) / 2
        T[0, 0, 1] = h * (kx ** (-1)) * sin(L * (kx)) - (h * (1 - cos(L * (kx))) * (kx ** (-1)) * sin(L * (kx))) / 3 + (
                    2 * h ** 3 * (1 - cos(L * (kx))) * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * sin(
                L * (kx))) / (3 * (kx ** 3))
        T[0, 0, 5] = (h ** 2 * L * (kx ** (-1)) * (1 + (kx ** 2) / h ** 2) * sin(L * (kx))) / 2 - h ** 2 * (
                    (1 - cos(L * (kx))) * (kx ** (-2)) + (
                        -((1 - cos(L * (kx))) * (kx ** (-2))) - (kx ** (-2)) * sin(L * (kx)) ** 2) / 3) + 2 * h ** 4 * (
                                 kx ** (-2)) * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 (L * (kx ** (-1)) * sin(L * (kx))) / 2 + (
                                     -((1 - cos(L * (kx))) * (kx ** (-2))) - (kx ** (-2)) * sin(L * (kx)) ** 2) / 3)
        T[0, 1, 1] = (h * ((1 - cos(L * (kx))) * (kx ** (-2)) + (kx ** (-2)) * sin(L * (kx)) ** 2)) / 6 + h ** 3 * (
                    kx ** (-2)) * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 (1 - cos(L * (kx))) * (kx ** (-2)) + (
                                     -((1 - cos(L * (kx))) * (kx ** (-2))) - (kx ** (-2)) * sin(L * (kx)) ** 2) / 3)
        T[0, 1, 5] = (h ** 2 * (1 - cos(L * (kx))) * sin(L * (kx))) / (3 * (kx ** 3)) + (
                    h ** 2 * (kx ** (-2)) * (1 + (kx ** 2) / h ** 2) * (
                        -(L * cos(L * (kx))) + (kx ** (-1)) * sin(L * (kx)))) / 2 + 2 * h ** 4 * (kx ** (-2)) * (
                                 -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 -((1 - cos(L * (kx))) * sin(L * (kx))) / (3 * (kx ** 3)) + (
                                     (kx ** (-2)) * (-(L * cos(L * (kx))) + (kx ** (-1)) * sin(L * (kx)))) / 2)
        T[0, 5, 5] = -(h * (1 - cos(L * (kx))) * (kx ** (-2))) + h ** 3 * (kx ** (-2)) * (1 + (kx ** 2) / h ** 2) * (
                    (1 - cos(L * (kx))) * (kx ** (-2)) - (L * (kx ** (-1)) * sin(L * (kx))) / 2) + (
                                 h ** 5 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                     (4 * (1 - cos(L * (kx))) * (kx ** (-2))) / 3 - L * (kx ** (-1)) * sin(L * (kx)) + (
                                         (kx ** (-2)) * sin(L * (kx)) ** 2) / 3)) / (kx ** 2) ** 2 + (
                                 h ** 3 * (kx ** (-2)) * ((1 - cos(L * (kx))) * (kx ** (-2)) + (
                                     -((1 - cos(L * (kx))) * (kx ** (-2))) - (kx ** (-2)) * sin(
                                 L * (kx)) ** 2) / 3)) / 2
        T[0, 2, 2] = (h * K1 * (1 - cos(L * (kx))) * (kx ** (-2))) / 2 + K2 * ((1 - cos(L * (kx))) * (kx ** (-2)) + (
                    K1 * (-2 * (1 - cos(L * (kx))) * (kx ** (-2)) - sin(sqrt(-K1) * L) ** 2 / K1)) / (h ** 2 + 5 * K1))
        T[0, 2, 3] = (2 * K2 * (
                    (cos(sqrt(-K1) * L) * sin(sqrt(-K1) * L)) / sqrt(-K1) - (kx ** (-1)) * sin(L * (kx)))) / (
                                 h ** 2 + 5 * K1)
        T[0, 3, 3] = -(h * (1 - cos(L * (kx))) * (kx ** (-2))) / 2 + (
                    K2 * (-2 * (1 - cos(L * (kx))) * (kx ** (-2)) - sin(sqrt(-K1) * L) ** 2 / K1)) / (h ** 2 + 5 * K1)
        T[1, 0, 0] = (h * (1 - cos(L * (kx))) * (kx) * sin(L * (kx))) / 3 + h * cos(L * (kx)) * (kx) * sin(L * (kx)) + (
                    h ** 3 * (1 + 2 * cos(L * (kx))) * (kx ** (-1)) * (
                        -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * sin(L * (kx))) / 3
        T[1, 0, 1] = h * cos(L * (kx)) - h * (cos(L * (kx)) ** 2 - sin(L * (kx)) ** 2) - (h * (kx ** 2) * (
                    -((1 - cos(L * (kx))) * (kx ** (-2))) + 2 * (kx ** (-2)) * sin(L * (kx)) ** 2)) / 3 + (
                                 2 * h ** 3 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                     -((1 - cos(L * (kx))) * (kx ** (-2))) + 2 * (kx ** (-2)) * sin(L * (kx)) ** 2)) / 3
        T[1, 0, 5] = (-2 * h ** 2 * (1 - cos(L * (kx))) * (kx ** (-1)) * sin(L * (kx))) / 3 + (
                    h ** 2 * (1 + (kx ** 2) / h ** 2) * (
                        L * cos(L * (kx)) + (kx ** (-1)) * sin(L * (kx)))) / 2 + 2 * h ** 4 * (kx ** (-2)) * (
                                 -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 (L * cos(L * (kx))) / 2 + ((kx ** (-1)) * sin(L * (kx))) / 6 - (
                                     2 * cos(L * (kx)) * (kx ** (-1)) * sin(L * (kx))) / 3) - h * (
                                 -(h * (1 - cos(L * (kx))) * (kx ** (-1)) * sin(L * (kx))) + h * cos(L * (kx)) * (
                                     kx ** (-1)) * sin(L * (kx)))
        T[1, 1, 1] = -(h * cos(L * (kx)) * (kx ** (-1)) * sin(L * (kx))) + (
                    h * (1 + 2 * cos(L * (kx))) * (kx ** (-1)) * sin(L * (kx))) / 6 + (
                                 2 * h ** 3 * (1 - cos(L * (kx))) * (
                                     -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * sin(L * (kx))) / (3 * (kx ** 3))
        T[1, 1, 5] = (h ** 2 * L * (kx ** (-1)) * (1 + (kx ** 2) / h ** 2) * sin(L * (kx))) / 2 + 2 * h ** 4 * (
                    kx ** (-2)) * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 ((1 - cos(L * (kx))) * (kx ** (-2))) / 3 + (L * (kx ** (-1)) * sin(L * (kx))) / 2 - (
                                     2 * (kx ** (-2)) * sin(L * (kx)) ** 2) / 3) + (h ** 2 * (
                    -((1 - cos(L * (kx))) * (kx ** (-2))) + 2 * (kx ** (-2)) * sin(L * (kx)) ** 2)) / 3 - h * (
                                 h * (1 - cos(L * (kx))) * cos(L * (kx)) * (kx ** (-2)) + h * (kx ** (-2)) * sin(
                             L * (kx)) ** 2)
        T[1, 5, 5] = -(h * (kx ** (-1)) * sin(L * (kx))) - (2 * h ** 3 * (1 - cos(L * (kx))) * sin(L * (kx))) / (
                    3 * (kx ** 3)) + (h ** 3 * (kx ** (-2)) * (1 + (kx ** 2) / h ** 2) * (
                    -(L * cos(L * (kx))) + (kx ** (-1)) * sin(L * (kx)))) / 2 + (
                                 h ** 5 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                     -(L * cos(L * (kx))) + ((kx ** (-1)) * sin(L * (kx))) / 3 + (
                                         2 * cos(L * (kx)) * (kx ** (-1)) * sin(L * (kx))) / 3)) / (kx ** 2) ** 2
        T[1, 2, 2] = (h * K1 * (kx ** (-1)) * sin(L * (kx))) / 2 + K2 * ((kx ** (-1)) * sin(L * (kx)) + (
                    2 * K1 * ((cos(sqrt(-K1) * L) * sin(sqrt(-K1) * L)) / sqrt(-K1) - (kx ** (-1)) * sin(L * (kx)))) / (
                                                                                     h ** 2 + 5 * K1))
        T[1, 2, 3] = (2 * K2 * (1 - cos(L * (kx)) - 2 * sin(sqrt(-K1) * L) ** 2)) / (h ** 2 + 5 * K1)
        T[1, 3, 3] = -(h * (kx ** (-1)) * sin(L * (kx))) / 2 + (
                    2 * K2 * ((cos(sqrt(-K1) * L) * sin(sqrt(-K1) * L)) / sqrt(-K1) - (kx ** (-1)) * sin(L * (kx)))) / (
                                 h ** 2 + 5 * K1)
        T[2, 0, 2] = -((h * K1 * (kx ** 2) * (
                    2 * cos(sqrt(-K1) * L) * (1 - cos(L * (kx))) * (kx ** (-2)) - ((sqrt(-K1) * kx) ** (-1)) * sin(
                sqrt(-K1) * L) * sin(L * (kx)))) / (h ** 2 + 5 * K1)) + (
                                 2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                                     cos(sqrt(-K1) * L) * (1 - cos(L * (kx))) + 2 * K1 * (
                                         (sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sin(L * (kx)))) / (
                                 h ** 2 + 5 * K1)
        T[2, 0, 3] = (h * sin(sqrt(-K1) * L)) / sqrt(-K1) + (2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                    -(((1 + cos(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + 2 * cos(sqrt(-K1) * L) * (
                        kx ** (-1)) * sin(L * (kx)))) / (h ** 2 + 5 * K1) - h * (kx ** 2) * (
                                 ((kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + (
                                     (-2 * K1 * (1 + cos(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(
                                 -K1) - cos(sqrt(-K1) * L) * (kx ** (-1)) * sin(L * (kx))) / (h ** 2 + 5 * K1))
        T[2, 1, 2] = (h * K1 * (-(((1 + cos(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + 2 * cos(sqrt(-K1) * L) * (
                    kx ** (-1)) * sin(L * (kx)))) / (h ** 2 + 5 * K1) + 2 * h ** 3 * (K1 / h ** 2 + K2 / h ** 3) * (
                                 ((kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + (
                                     (-2 * K1 * (1 + cos(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(
                                 -K1) - cos(sqrt(-K1) * L) * (kx ** (-1)) * sin(L * (kx))) / (h ** 2 + 5 * K1))
        T[2, 1, 3] = (2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                    2 * cos(sqrt(-K1) * L) * (1 - cos(L * (kx))) * (kx ** (-2)) - ((sqrt(-K1) * kx) ** (-1)) * sin(
                sqrt(-K1) * L) * sin(L * (kx)))) / (h ** 2 + 5 * K1) + (h * (
                    cos(sqrt(-K1) * L) * (1 - cos(L * (kx))) + 2 * K1 * ((sqrt(-K1) * kx) ** (-1)) * sin(
                sqrt(-K1) * L) * sin(L * (kx)))) / (h ** 2 + 5 * K1)
        T[2, 2, 5] = -(K1 * L * sin(sqrt(-K1) * L)) / (2 * sqrt(-K1)) + (h ** 2 * K1 * (
                    2 * cos(sqrt(-K1) * L) * (1 - cos(L * (kx))) * (kx ** (-2)) - ((sqrt(-K1) * kx) ** (-1)) * sin(
                sqrt(-K1) * L) * sin(L * (kx)))) / (h ** 2 + 5 * K1) + 2 * h ** 4 * (kx ** (-2)) * (
                                 -1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                                 (L * sin(sqrt(-K1) * L)) / (2 * sqrt(-K1)) - (
                                     cos(sqrt(-K1) * L) * (1 - cos(L * (kx))) + 2 * K1 * (
                                         (sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sin(L * (kx))) / (
                                             h ** 2 + 5 * K1))
        T[2, 3, 5] = (-(L * cos(sqrt(-K1) * L)) + sin(sqrt(-K1) * L) / sqrt(-K1)) / 2 + h ** 2 * (
                    ((kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + (
                        (-2 * K1 * (1 + cos(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) - cos(
                    sqrt(-K1) * L) * (kx ** (-1)) * sin(L * (kx))) / (h ** 2 + 5 * K1)) + 2 * h ** 4 * (kx ** (-2)) * (
                                 -1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                                 -(-(L * cos(sqrt(-K1) * L)) + sin(sqrt(-K1) * L) / sqrt(-K1)) / (2 * K1) - (
                                     -(((1 + cos(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + 2 * cos(
                                 sqrt(-K1) * L) * (kx ** (-1)) * sin(L * (kx))) / (h ** 2 + 5 * K1))
        T[3, 0, 2] = -((h * K1 * cos(L * (kx)) * sin(sqrt(-K1) * L)) / sqrt(-K1)) - (h * K1 * (kx ** 2) * (
                    -((cos(L * (kx)) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + (
                        2 * K1 * (1 - cos(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + cos(
                sqrt(-K1) * L) * (kx ** (-1)) * sin(L * (kx)))) / (h ** 2 + 5 * K1) + (
                                 2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                                     (K1 * (1 + cos(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1) + (
                                         h ** 2 + 3 * K1) * cos(sqrt(-K1) * L) * (kx ** (-1)) * sin(L * (kx)))) / (
                                 h ** 2 + 5 * K1)
        T[3, 0, 3] = h * cos(sqrt(-K1) * L) - h * cos(sqrt(-K1) * L) * cos(L * (kx)) + (
                    2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                        -(cos(sqrt(-K1) * L) * (1 - cos(L * (kx)))) + (h ** 2 + 3 * K1) * (
                            (sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sin(L * (kx)))) / (h ** 2 + 5 * K1) - h * (
                                 kx ** 2) * (cos(sqrt(-K1) * L) * (kx ** (-2)) + (
                    -(cos(sqrt(-K1) * L) * cos(L * (kx))) - 2 * K1 * cos(sqrt(-K1) * L) * (1 + cos(L * (kx))) * (
                        kx ** (-2)) + K1 * ((sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sin(L * (kx))) / (
                                                         h ** 2 + 5 * K1))
        T[3, 1, 2] = -(h * K1 * ((sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sin(L * (kx))) + (h * K1 * (
                    -(cos(sqrt(-K1) * L) * (1 - cos(L * (kx)))) + (h ** 2 + 3 * K1) * ((sqrt(-K1) * kx) ** (-1)) * sin(
                sqrt(-K1) * L) * sin(L * (kx)))) / (h ** 2 + 5 * K1) + 2 * h ** 3 * (K1 / h ** 2 + K2 / h ** 3) * (
                                 cos(sqrt(-K1) * L) * (kx ** (-2)) + (
                                     -(cos(sqrt(-K1) * L) * cos(L * (kx))) - 2 * K1 * cos(sqrt(-K1) * L) * (
                                         1 + cos(L * (kx))) * (kx ** (-2)) + K1 * ((sqrt(-K1) * kx) ** (-1)) * sin(
                                 sqrt(-K1) * L) * sin(L * (kx))) / (h ** 2 + 5 * K1))
        T[3, 1, 3] = -(h * cos(sqrt(-K1) * L) * (kx ** (-1)) * sin(L * (kx))) + (
                    2 * h ** 3 * (K1 / h ** 2 + K2 / h ** 3) * (-((cos(L * (kx)) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + (
                        2 * K1 * (1 - cos(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + cos(
                sqrt(-K1) * L) * (kx ** (-1)) * sin(L * (kx)))) / (h ** 2 + 5 * K1) + (h * (
                    (K1 * (1 + cos(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1) + (h ** 2 + 3 * K1) * cos(
                sqrt(-K1) * L) * (kx ** (-1)) * sin(L * (kx)))) / (h ** 2 + 5 * K1)
        T[3, 2, 5] = -((h ** 2 * K1 * (1 - cos(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1)) - (
                    K1 * (L * cos(sqrt(-K1) * L) + sin(sqrt(-K1) * L) / sqrt(-K1))) / 2 + (h ** 2 * K1 * (
                    -((cos(L * (kx)) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + (
                        2 * K1 * (1 - cos(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + cos(
                sqrt(-K1) * L) * (kx ** (-1)) * sin(L * (kx)))) / (h ** 2 + 5 * K1) + 2 * h ** 4 * (
                                 K1 / h ** 2 + K2 / h ** 3) * (kx ** (-2)) * (
                                 (L * cos(sqrt(-K1) * L)) / 2 + sin(sqrt(-K1) * L) / (2 * sqrt(-K1)) + (
                                     -((K1 * (1 + cos(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1)) - (
                                         h ** 2 + 3 * K1) * cos(sqrt(-K1) * L) * (kx ** (-1)) * sin(L * (kx))) / (
                                             h ** 2 + 5 * K1))
        T[3, 3, 5] = -(h ** 2 * cos(sqrt(-K1) * L) * (1 - cos(L * (kx))) * (kx ** (-2))) - (
                    K1 * L * sin(sqrt(-K1) * L)) / (2 * sqrt(-K1)) + h ** 2 * (cos(sqrt(-K1) * L) * (kx ** (-2)) + (
                    -(cos(sqrt(-K1) * L) * cos(L * (kx))) - 2 * K1 * cos(sqrt(-K1) * L) * (1 + cos(L * (kx))) * (
                        kx ** (-2)) + K1 * ((sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sin(L * (kx))) / (
                                                                                           h ** 2 + 5 * K1)) + 2 * h ** 4 * (
                                 K1 / h ** 2 + K2 / h ** 3) * (kx ** (-2)) * (
                                 (L * sin(sqrt(-K1) * L)) / (2 * sqrt(-K1)) - (
                                     -(cos(sqrt(-K1) * L) * (1 - cos(L * (kx)))) + (h ** 2 + 3 * K1) * (
                                         (sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sin(L * (kx))) / (
                                             h ** 2 + 5 * K1))
        return T

    if h ** 2 + K1 < 0.0 and h > 1e-6:
        kx = sqrt(-kx2)
        ky = sqrt(ky2)

        T[0, 0, 0] = (h ** 3 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                    (1 - cosh(L * (kx))) * (kx ** (-2)) - (kx ** (-2)) * sinh(L * (kx)) ** 2)) / 3 + (h * (kx ** 2) * (
                    (1 - cosh(L * (kx))) * (kx ** (-2)) + (
                        -((1 - cosh(L * (kx))) * (kx ** (-2))) + (kx ** (-2)) * sinh(L * (kx)) ** 2) / 3)) / 2
        T[0, 0, 1] = h * (kx ** (-1)) * sinh(L * (kx)) - (
                    h * (1 - cosh(L * (kx))) * (kx ** (-1)) * sinh(L * (kx))) / 3 + (
                                 2 * h ** 3 * (1 - cosh(L * (kx))) * (
                                     -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * sinh(L * (kx))) / (
                                 3 * (kx ** 3))
        T[0, 0, 5] = (h ** 2 * L * (kx ** (-1)) * (1 + (kx ** 2) / h ** 2) * sinh(L * (kx))) / 2 - h ** 2 * (
                    (1 - cosh(L * (kx))) * (kx ** (-2)) + (-((1 - cosh(L * (kx))) * (kx ** (-2))) + (kx ** (-2)) * sinh(
                L * (kx)) ** 2) / 3) + 2 * h ** 4 * (kx ** (-2)) * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 (L * (kx ** (-1)) * sinh(L * (kx))) / 2 + (
                                     -((1 - cosh(L * (kx))) * (kx ** (-2))) + (kx ** (-2)) * sinh(L * (kx)) ** 2) / 3)
        T[0, 1, 1] = (h * ((1 - cosh(L * (kx))) * (kx ** (-2)) - (kx ** (-2)) * sinh(L * (kx)) ** 2)) / 6 + h ** 3 * (
                    kx ** (-2)) * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 (1 - cosh(L * (kx))) * (kx ** (-2)) + (
                                     -((1 - cosh(L * (kx))) * (kx ** (-2))) + (kx ** (-2)) * sinh(L * (kx)) ** 2) / 3)
        T[0, 1, 5] = (h ** 2 * (1 - cosh(L * (kx))) * sinh(L * (kx))) / (3 * (kx ** 3)) + (
                    h ** 2 * (kx ** (-2)) * (1 + (kx ** 2) / h ** 2) * (
                        -(L * cosh(L * (kx))) + (kx ** (-1)) * sinh(L * (kx)))) / 2 + 2 * h ** 4 * (kx ** (-2)) * (
                                 -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 -((1 - cosh(L * (kx))) * sinh(L * (kx))) / (3 * (kx ** 3)) + (
                                     (kx ** (-2)) * (-(L * cosh(L * (kx))) + (kx ** (-1)) * sinh(L * (kx)))) / 2)
        T[0, 5, 5] = -(h * (1 - cosh(L * (kx))) * (kx ** (-2))) + h ** 3 * (kx ** (-2)) * (1 + (kx ** 2) / h ** 2) * (
                    (1 - cosh(L * (kx))) * (kx ** (-2)) - (L * (kx ** (-1)) * sinh(L * (kx))) / 2) + (
                                 h ** 5 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                     (4 * (1 - cosh(L * (kx))) * (kx ** (-2))) / 3 - L * (kx ** (-1)) * sinh(
                                 L * (kx)) - ((kx ** (-2)) * sinh(L * (kx)) ** 2) / 3)) / (kx ** 2) ** 2 + (
                                 h ** 3 * (kx ** (-2)) * ((1 - cosh(L * (kx))) * (kx ** (-2)) + (
                                     -((1 - cosh(L * (kx))) * (kx ** (-2))) + (kx ** (-2)) * sinh(
                                 L * (kx)) ** 2) / 3)) / 2
        T[0, 2, 2] = (h * K1 * (1 - cosh(L * (kx))) * (kx ** (-2))) / 2 + K2 * ((1 - cosh(L * (kx))) * (kx ** (-2)) + (
                    K1 * (-2 * (1 - cosh(L * (kx))) * (kx ** (-2)) - sin(sqrt(-K1) * L) ** 2 / K1)) / (h ** 2 + 5 * K1))
        T[0, 2, 3] = (2 * K2 * (
                    (cos(sqrt(-K1) * L) * sin(sqrt(-K1) * L)) / sqrt(-K1) - (kx ** (-1)) * sinh(L * (kx)))) / (
                                 h ** 2 + 5 * K1)
        T[0, 3, 3] = -(h * (1 - cosh(L * (kx))) * (kx ** (-2))) / 2 + (
                    K2 * (-2 * (1 - cosh(L * (kx))) * (kx ** (-2)) - sin(sqrt(-K1) * L) ** 2 / K1)) / (h ** 2 + 5 * K1)
        T[1, 0, 0] = -(h * (1 - cosh(L * (kx))) * (kx) * sinh(L * (kx))) / 3 - h * cosh(L * (kx)) * (kx) * sinh(
            L * (kx)) + (h ** 3 * (1 + 2 * cosh(L * (kx))) * (kx ** (-1)) * (
                    -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * sinh(L * (kx))) / 3
        T[1, 0, 1] = h * cosh(L * (kx)) - h * (cosh(L * (kx)) ** 2 + sinh(L * (kx)) ** 2) - (h * (kx ** 2) * (
                    -((1 - cosh(L * (kx))) * (kx ** (-2))) - 2 * (kx ** (-2)) * sinh(L * (kx)) ** 2)) / 3 + (
                                 2 * h ** 3 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                     -((1 - cosh(L * (kx))) * (kx ** (-2))) - 2 * (kx ** (-2)) * sinh(
                                 L * (kx)) ** 2)) / 3
        T[1, 0, 5] = (-2 * h ** 2 * (1 - cosh(L * (kx))) * (kx ** (-1)) * sinh(L * (kx))) / 3 + (
                    h ** 2 * (1 + (kx ** 2) / h ** 2) * (
                        L * cosh(L * (kx)) + (kx ** (-1)) * sinh(L * (kx)))) / 2 + 2 * h ** 4 * (kx ** (-2)) * (
                                 -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 (L * cosh(L * (kx))) / 2 + ((kx ** (-1)) * sinh(L * (kx))) / 6 - (
                                     2 * cosh(L * (kx)) * (kx ** (-1)) * sinh(L * (kx))) / 3) - h * (
                                 -(h * (1 - cosh(L * (kx))) * (kx ** (-1)) * sinh(L * (kx))) + h * cosh(L * (kx)) * (
                                     kx ** (-1)) * sinh(L * (kx)))
        T[1, 1, 1] = -(h * cosh(L * (kx)) * (kx ** (-1)) * sinh(L * (kx))) + (
                    h * (1 + 2 * cosh(L * (kx))) * (kx ** (-1)) * sinh(L * (kx))) / 6 + (
                                 2 * h ** 3 * (1 - cosh(L * (kx))) * (
                                     -1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * sinh(L * (kx))) / (
                                 3 * (kx ** 3))
        T[1, 1, 5] = (h ** 2 * L * (kx ** (-1)) * (1 + (kx ** 2) / h ** 2) * sinh(L * (kx))) / 2 + (h ** 2 * (
                    -((1 - cosh(L * (kx))) * (kx ** (-2))) - 2 * (kx ** (-2)) * sinh(
                L * (kx)) ** 2)) / 3 + 2 * h ** 4 * (kx ** (-2)) * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                 ((1 - cosh(L * (kx))) * (kx ** (-2))) / 3 + (L * (kx ** (-1)) * sinh(L * (kx))) / 2 + (
                                     2 * (kx ** (-2)) * sinh(L * (kx)) ** 2) / 3) - h * (
                                 h * (1 - cosh(L * (kx))) * cosh(L * (kx)) * (kx ** (-2)) - h * (kx ** (-2)) * sinh(
                             L * (kx)) ** 2)
        T[1, 5, 5] = -(h * (kx ** (-1)) * sinh(L * (kx))) - (2 * h ** 3 * (1 - cosh(L * (kx))) * sinh(L * (kx))) / (
                    3 * (kx ** 3)) + (h ** 3 * (kx ** (-2)) * (1 + (kx ** 2) / h ** 2) * (
                    -(L * cosh(L * (kx))) + (kx ** (-1)) * sinh(L * (kx)))) / 2 + (
                                 h ** 5 * (-1 - K2 / h ** 3 + 2 * (1 - (kx ** 2) / h ** 2)) * (
                                     -(L * cosh(L * (kx))) + ((kx ** (-1)) * sinh(L * (kx))) / 3 + (
                                         2 * cosh(L * (kx)) * (kx ** (-1)) * sinh(L * (kx))) / 3)) / (kx ** 2) ** 2
        T[1, 2, 2] = (h * K1 * (kx ** (-1)) * sinh(L * (kx))) / 2 + K2 * ((kx ** (-1)) * sinh(L * (kx)) + (2 * K1 * (
                    (cos(sqrt(-K1) * L) * sin(sqrt(-K1) * L)) / sqrt(-K1) - (kx ** (-1)) * sinh(L * (kx)))) / (
                                                                                      h ** 2 + 5 * K1))
        T[1, 2, 3] = (2 * K2 * (1 - cosh(L * (kx)) - 2 * sin(sqrt(-K1) * L) ** 2)) / (h ** 2 + 5 * K1)
        T[1, 3, 3] = -(h * (kx ** (-1)) * sinh(L * (kx))) / 2 + (2 * K2 * (
                    (cos(sqrt(-K1) * L) * sin(sqrt(-K1) * L)) / sqrt(-K1) - (kx ** (-1)) * sinh(L * (kx)))) / (
                                 h ** 2 + 5 * K1)
        T[2, 0, 2] = -((h * K1 * (kx ** 2) * (
                    2 * cos(sqrt(-K1) * L) * (1 - cosh(L * (kx))) * (kx ** (-2)) - ((sqrt(-K1) * kx) ** (-1)) * sin(
                sqrt(-K1) * L) * sinh(L * (kx)))) / (h ** 2 + 5 * K1)) + (
                                 2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                                     cos(sqrt(-K1) * L) * (1 - cosh(L * (kx))) + 2 * K1 * (
                                         (sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sinh(L * (kx)))) / (
                                 h ** 2 + 5 * K1)
        T[2, 0, 3] = (h * sin(sqrt(-K1) * L)) / sqrt(-K1) + (2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                    -(((1 + cosh(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + 2 * cos(sqrt(-K1) * L) * (
                        kx ** (-1)) * sinh(L * (kx)))) / (h ** 2 + 5 * K1) - h * (kx ** 2) * (
                                 ((kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + (
                                     (-2 * K1 * (1 + cosh(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(
                                 -K1) - cos(sqrt(-K1) * L) * (kx ** (-1)) * sinh(L * (kx))) / (h ** 2 + 5 * K1))
        T[2, 1, 2] = (h * K1 * (-(((1 + cosh(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + 2 * cos(sqrt(-K1) * L) * (
                    kx ** (-1)) * sinh(L * (kx)))) / (h ** 2 + 5 * K1) + 2 * h ** 3 * (K1 / h ** 2 + K2 / h ** 3) * (
                                 ((kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + (
                                     (-2 * K1 * (1 + cosh(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(
                                 -K1) - cos(sqrt(-K1) * L) * (kx ** (-1)) * sinh(L * (kx))) / (h ** 2 + 5 * K1))
        T[2, 1, 3] = (2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                    2 * cos(sqrt(-K1) * L) * (1 - cosh(L * (kx))) * (kx ** (-2)) - ((sqrt(-K1) * kx) ** (-1)) * sin(
                sqrt(-K1) * L) * sinh(L * (kx)))) / (h ** 2 + 5 * K1) + (h * (
                    cos(sqrt(-K1) * L) * (1 - cosh(L * (kx))) + 2 * K1 * ((sqrt(-K1) * kx) ** (-1)) * sin(
                sqrt(-K1) * L) * sinh(L * (kx)))) / (h ** 2 + 5 * K1)
        T[2, 2, 5] = -(K1 * L * sin(sqrt(-K1) * L)) / (2 * sqrt(-K1)) + (h ** 2 * K1 * (
                    2 * cos(sqrt(-K1) * L) * (1 - cosh(L * (kx))) * (kx ** (-2)) - ((sqrt(-K1) * kx) ** (-1)) * sin(
                sqrt(-K1) * L) * sinh(L * (kx)))) / (h ** 2 + 5 * K1) + 2 * h ** 4 * (kx ** (-2)) * (
                                 -1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                                 (L * sin(sqrt(-K1) * L)) / (2 * sqrt(-K1)) - (
                                     cos(sqrt(-K1) * L) * (1 - cosh(L * (kx))) + 2 * K1 * (
                                         (sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sinh(L * (kx))) / (
                                             h ** 2 + 5 * K1))
        T[2, 3, 5] = (-(L * cos(sqrt(-K1) * L)) + sin(sqrt(-K1) * L) / sqrt(-K1)) / 2 + h ** 2 * (
                    ((kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + (
                        (-2 * K1 * (1 + cosh(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) - cos(
                    sqrt(-K1) * L) * (kx ** (-1)) * sinh(L * (kx))) / (h ** 2 + 5 * K1)) + 2 * h ** 4 * (kx ** (-2)) * (
                                 -1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                                 -(-(L * cos(sqrt(-K1) * L)) + sin(sqrt(-K1) * L) / sqrt(-K1)) / (2 * K1) - (
                                     -(((1 + cosh(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + 2 * cos(
                                 sqrt(-K1) * L) * (kx ** (-1)) * sinh(L * (kx))) / (h ** 2 + 5 * K1))
        T[3, 0, 2] = -((h * K1 * cosh(L * (kx)) * sin(sqrt(-K1) * L)) / sqrt(-K1)) - (h * K1 * (kx ** 2) * (
                    -((cosh(L * (kx)) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + (
                        2 * K1 * (1 - cosh(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + cos(
                sqrt(-K1) * L) * (kx ** (-1)) * sinh(L * (kx)))) / (h ** 2 + 5 * K1) + (
                                 2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                                     (K1 * (1 + cosh(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1) + (
                                         h ** 2 + 3 * K1) * cos(sqrt(-K1) * L) * (kx ** (-1)) * sinh(L * (kx)))) / (
                                 h ** 2 + 5 * K1)
        T[3, 0, 3] = h * cos(sqrt(-K1) * L) - h * cos(sqrt(-K1) * L) * cosh(L * (kx)) + (
                    2 * h ** 3 * (-1 + K2 / h ** 3 + (kx ** 2) / h ** 2) * (
                        -(cos(sqrt(-K1) * L) * (1 - cosh(L * (kx)))) + (h ** 2 + 3 * K1) * (
                            (sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sinh(L * (kx)))) / (
                                 h ** 2 + 5 * K1) - h * (kx ** 2) * (cos(sqrt(-K1) * L) * (kx ** (-2)) + (
                    -(cos(sqrt(-K1) * L) * cosh(L * (kx))) - 2 * K1 * cos(sqrt(-K1) * L) * (1 + cosh(L * (kx))) * (
                        kx ** (-2)) + K1 * ((sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sinh(L * (kx))) / (
                                                                                 h ** 2 + 5 * K1))
        T[3, 1, 2] = -(h * K1 * ((sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sinh(L * (kx))) + (h * K1 * (
                    -(cos(sqrt(-K1) * L) * (1 - cosh(L * (kx)))) + (h ** 2 + 3 * K1) * ((sqrt(-K1) * kx) ** (-1)) * sin(
                sqrt(-K1) * L) * sinh(L * (kx)))) / (h ** 2 + 5 * K1) + 2 * h ** 3 * (K1 / h ** 2 + K2 / h ** 3) * (
                                 cos(sqrt(-K1) * L) * (kx ** (-2)) + (
                                     -(cos(sqrt(-K1) * L) * cosh(L * (kx))) - 2 * K1 * cos(sqrt(-K1) * L) * (
                                         1 + cosh(L * (kx))) * (kx ** (-2)) + K1 * ((sqrt(-K1) * kx) ** (-1)) * sin(
                                 sqrt(-K1) * L) * sinh(L * (kx))) / (h ** 2 + 5 * K1))
        T[3, 1, 3] = -(h * cos(sqrt(-K1) * L) * (kx ** (-1)) * sinh(L * (kx))) + (
                    2 * h ** 3 * (K1 / h ** 2 + K2 / h ** 3) * (-((cosh(L * (kx)) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + (
                        2 * K1 * (1 - cosh(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + cos(
                sqrt(-K1) * L) * (kx ** (-1)) * sinh(L * (kx)))) / (h ** 2 + 5 * K1) + (h * (
                    (K1 * (1 + cosh(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1) + (h ** 2 + 3 * K1) * cos(
                sqrt(-K1) * L) * (kx ** (-1)) * sinh(L * (kx)))) / (h ** 2 + 5 * K1)
        T[3, 2, 5] = -((h ** 2 * K1 * (1 - cosh(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1)) - (
                    K1 * (L * cos(sqrt(-K1) * L) + sin(sqrt(-K1) * L) / sqrt(-K1))) / 2 + (h ** 2 * K1 * (
                    -
((cosh(L * (kx)) * sin(sqrt(-K1) * L)) / sqrt(-K1)) + (
                        2 * K1 * (1 - cosh(L * (kx))) * (kx ** (-2)) * sin(sqrt(-K1) * L)) / sqrt(-K1) + cos(
                sqrt(-K1) * L) * (kx ** (-1)) * sinh(L * (kx)))) / (h ** 2 + 5 * K1) + 2 * h ** 4 * (
                                 K1 / h ** 2 + K2 / h ** 3) * (kx ** (-2)) * (
                                 (L * cos(sqrt(-K1) * L)) / 2 + sin(sqrt(-K1) * L) / (2 * sqrt(-K1)) + (
                                     -((K1 * (1 + cosh(L * (kx))) * sin(sqrt(-K1) * L)) / sqrt(-K1)) - (
                                         h ** 2 + 3 * K1) * cos(sqrt(-K1) * L) * (kx ** (-1)) * sinh(L * (kx))) / (
                                             h ** 2 + 5 * K1))
        T[3, 3, 5] = -(h ** 2 * cos(sqrt(-K1) * L) * (1 - cosh(L * (kx))) * (kx ** (-2))) - (
                    K1 * L * sin(sqrt(-K1) * L)) / (2 * sqrt(-K1)) + h ** 2 * (cos(sqrt(-K1) * L) * (kx ** (-2)) + (
                    -(cos(sqrt(-K1) * L) * cosh(L * (kx))) - 2 * K1 * cos(sqrt(-K1) * L) * (1 + cosh(L * (kx))) * (
                        kx ** (-2)) + K1 * ((sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sinh(L * (kx))) / (
                                                                                           h ** 2 + 5 * K1)) + 2 * h ** 4 * (
                                 K1 / h ** 2 + K2 / h ** 3) * (kx ** (-2)) * (
                                 (L * sin(sqrt(-K1) * L)) / (2 * sqrt(-K1)) - (
                                     -(cos(sqrt(-K1) * L) * (1 - cosh(L * (kx)))) + (h ** 2 + 3 * K1) * (
                                         (sqrt(-K1) * kx) ** (-1)) * sin(sqrt(-K1) * L) * sinh(L * (kx))) / (
                                             h ** 2 + 5 * K1))
        return T

    if np.abs(kx2 - 4 * ky2) < 1e-3 and K1 != 0.0:

        T[0, 0, 0] = (-3 * h ** 3 * L ** 2) / 10 - (K2 * L ** 2) / 2 + (13 * h ** 5 * L ** 4) / 150 + (
                    h ** 2 * K2 * L ** 4) / 10 - (53 * h ** 7 * L ** 6) / 5625 - (11 * h ** 4 * K2 * L ** 6) / 1125 + (
                                 71 * h ** 9 * L ** 8) / 131250 + (43 * h ** 6 * K2 * L ** 8) / 78750 - (
                                 19 * h ** 8 * K2 * L ** 10) / 984375 + 4.6723317834428953 * 10 ** -7 * h ** 10 * K2 * L ** 12
        T[0, 0, 1] = h * L - (7 * h ** 3 * L ** 3) / 15 - (K2 * L ** 3) / 3 + (9 * h ** 5 * L ** 5) / 125 + (
                    h ** 2 * K2 * L ** 5) / 15 - (214 * h ** 7 * L ** 7) / 39375 - (2 * h ** 4 * K2 * L ** 7) / 375 + (
                                 61 * h ** 9 * L ** 9) / 253125 + (17 * h ** 6 * K2 * L ** 9) / 70875 - (
                                 62 * h ** 8 * K2 * L ** 11) / 8859375 + (2 * h ** 10 * K2 * L ** 13) / 13921875
        T[0, 0, 5] = (9 * h ** 2 * L ** 2) / 10 - (h * K2 * L ** 4) / 12 + (4 * h ** 3 * K2 * L ** 6) / 225 - (
                    13 * h ** 5 * K2 * L ** 8) / 10500 + (83 * h ** 7 * K2 * L ** 10) / 1771875 - (
                                 677 * h ** 9 * K2 * L ** 12) / 584718750 + (h ** 4 * (
                    (-18 * L ** 4) / 5 + (L * ((-9 * L ** 3) / sqrt(5) - 23 * sqrt(5) * L ** 3)) / sqrt(5))) / 120 + (
                                 h ** 6 * ((262 * L ** 6) / 125 - (
                                     L ** 3 * ((-9 * L ** 3) / sqrt(5) - 23 * sqrt(5) * L ** 3)) / (
                                                       30 * sqrt(5)))) / 120 + (h ** 8 * (
                    (-4091 * L ** 8) / 26250 + (L ** 5 * ((-9 * L ** 3) / sqrt(5) - 23 * sqrt(5) * L ** 3)) / (
                        3000 * sqrt(5)))) / 120 + (h ** 10 * (
                    (13603 * L ** 10) / 2362500 - (L ** 7 * ((-9 * L ** 3) / sqrt(5) - 23 * sqrt(5) * L ** 3)) / (
                        630000 * sqrt(5)))) / 120
        T[0, 1, 1] = (h * L ** 2) / 4 - (h ** 3 * L ** 4) / 10 - (K2 * L ** 4) / 12 + (13 * h ** 5 * L ** 6) / 1125 + (
                    h ** 2 * K2 * L ** 6) / 90 - (53 * h ** 7 * L ** 8) / 78750 - (h ** 4 * K2 * L ** 8) / 1500 + (
                                 71 * h ** 9 * L ** 10) / 2953125 + (17 * h ** 6 * K2 * L ** 10) / 708750 - (
                                 31 * h ** 8 * K2 * L ** 12) / 53156250 + (h ** 10 * K2 * L ** 14) / 97453125
        T[0, 1, 5] = (7 * h ** 2 * L ** 3) / 15 - (131 * h ** 4 * L ** 5) / 1500 - (h * K2 * L ** 5) / 20 + (
                    89 * h ** 6 * L ** 7) / 13125 + (h ** 3 * K2 * L ** 7) / 175 - (
                                 2137 * h ** 8 * L ** 9) / 7087500 - (h ** 5 * K2 * L ** 9) / 3500 + (
                                 122 * h ** 10 * L ** 11) / 13921875 + (8 * h ** 7 * K2 * L ** 11) / 928125 - (
                                 151 * h ** 9 * K2 * L ** 13) / 844593750
        T[0, 5, 5] = -(h * L ** 2) / 2 + (3 * h ** 3 * L ** 4) / 20 + 2.8912057932946805 * 10 ** -17 * K2 * L ** 4 - (
                    139 * h ** 5 * L ** 6) / 9000 - (h ** 2 * K2 * L ** 6) / 120 + (271 * h ** 7 * L ** 8) / 315000 + (
                                 h ** 4 * K2 * L ** 8) / 1400 - (143 * h ** 9 * L ** 10) / 4725000 - (
                                 h ** 6 * K2 * L ** 10) / 35000 + (
                                 2 * h ** 8 * K2 * L ** 12) / 2784375 - 1.27702984845842 * 10 ** -8 * h ** 10 * K2 * L ** 14
        T[0, 2, 2] = -(h ** 3 * L ** 2) / 20 + (K2 * L ** 2) / 2 + (h ** 5 * L ** 4) / 300 - (
                    h ** 2 * K2 * L ** 4) / 20 - (h ** 7 * L ** 6) / 11250 + (2 * h ** 4 * K2 * L ** 6) / 1125 + (
                                 h ** 9 * L ** 8) / 787500 - (h ** 6 * K2 * L ** 8) / 31500 + (
                                 h ** 8 * K2 * L ** 10) / 2953125 - (h ** 10 * K2 * L ** 12) / 417656250
        T[0, 2, 3] = (K2 * L ** 3) / 3 - (2 * h ** 2 * K2 * L ** 5) / 75 + (2 * h ** 4 * K2 * L ** 7) / 2625 - (
                    4 * h ** 6 * K2 * L ** 9) / 354375 + (2 * h ** 8 * K2 * L ** 11) / 19490625 - (
                                 4 * h ** 10 * K2 * L ** 13) / 6334453125
        T[0, 3, 3] = -(h * L ** 2) / 4 + (h ** 3 * L ** 4) / 60 + (K2 * L ** 4) / 12 - (h ** 5 * L ** 6) / 2250 - (
                    h ** 2 * K2 * L ** 6) / 225 + (h ** 7 * L ** 8) / 157500 + (h ** 4 * K2 * L ** 8) / 10500 - (
                                 h ** 9 * L ** 10) / 17718750 - (2 * h ** 6 * K2 * L ** 10) / 1771875 + (
                                 h ** 8 * K2 * L ** 12) / 116943750 - 4.5104807009568884 * 10 ** -11 * h ** 10 * K2 * L ** 14
        T[1, 0, 0] = (h ** 3 * L) / 5 - K2 * L - (2 * h ** 5 * L ** 3) / 25 + (2 * h ** 2 * K2 * L ** 3) / 5 + (
                    22 * h ** 7 * L ** 5) / 1875 - (22 * h ** 4 * K2 * L ** 5) / 375 - (
                                 172 * h ** 9 * L ** 7) / 196875 + (172 * h ** 6 * K2 * L ** 7) / 39375 - (
                                 38 * h ** 8 * K2 * L ** 9) / 196875 + 5.606798140131473 * 10 ** -6 * h ** 10 * K2 * L ** 11
        T[1, 0, 1] = (h ** 3 * L ** 2) / 5 - K2 * L ** 2 - (h ** 5 * L ** 4) / 15 + (h ** 2 * K2 * L ** 4) / 3 + (
                    14 * h ** 7 * L ** 6) / 1875 - (14 * h ** 4 * K2 * L ** 6) / 375 - (
                                 17 * h ** 9 * L ** 8) / 39375 + (17 * h ** 6 * K2 * L ** 8) / 7875 - (
                                 682 * h ** 8 * K2 * L ** 10) / 8859375 + (26 * h ** 10 * K2 * L ** 12) / 13921875
        T[1, 0, 5] = (4 * h ** 2 * L) / 5 - (h ** 4 * L ** 3) / 75 - (h * K2 * L ** 3) / 3 - (
                    12 * h ** 6 * L ** 5) / 625 + (8 * h ** 3 * K2 * L ** 5) / 75 + (386 * h ** 8 * L ** 7) / 196875 - (
                                 26 * h ** 5 * K2 * L ** 7) / 2625 - (166 * h ** 10 * L ** 9) / 1771875 + (
                                 166 * h ** 7 * K2 * L ** 9) / 354375 - (1354 * h ** 9 * K2 * L ** 11) / 97453125
        T[1, 1, 1] = -(h * L) / 2 + (2 * h ** 3 * L ** 3) / 15 - (K2 * L ** 3) / 3 - (2 * h ** 5 * L ** 5) / 125 + (
                    h ** 2 * K2 * L ** 5) / 15 + (44 * h ** 7 * L ** 7) / 39375 - (2 * h ** 4 * K2 * L ** 7) / 375 - (
                                 86 * h ** 9 * L ** 9) / 1771875 + (17 * h ** 6 * K2 * L ** 9) / 70875 - (
                                 62 * h ** 8 * K2 * L ** 11) / 8859375 + (2 * h ** 10 * K2 * L ** 13) / 13921875
        T[1, 1, 5] = -(h ** 2 * L ** 2) / 10 - (h * K2 * L ** 4) / 4 + (h ** 3 * K2 * L ** 6) / 25 - (
                    9 * h ** 5 * K2 * L ** 8) / 3500 + (8 * h ** 7 * K2 * L ** 10) / 84375 - (
                                 151 * h ** 9 * K2 * L ** 12) / 64968750 + (h ** 4 * (
                    L ** 4 / 5 - (L * ((9 * L ** 3) / (2 * sqrt(5)) - (9 * sqrt(5) * L ** 3) / 2)) / sqrt(5))) / 60 + (
                                 h ** 6 * ((-49 * L ** 6) / 125 + (
                                     L ** 3 * ((9 * L ** 3) / (2 * sqrt(5)) - (9 * sqrt(5) * L ** 3) / 2)) / (
                                                       30 * sqrt(5)))) / 60 + (h ** 8 * (
                    (227 * L ** 8) / 7500 - (L ** 5 * ((9 * L ** 3) / (2 * sqrt(5)) - (9 * sqrt(5) * L ** 3) / 2)) / (
                        3000 * sqrt(5)))) / 60 + (h ** 10 * ((-5381 * L ** 10) / 4725000 + (
                    L ** 7 * ((9 * L ** 3) / (2 * sqrt(5)) - (9 * sqrt(5) * L ** 3) / 2)) / (630000 * sqrt(5)))) / 60
        T[1, 5, 5] = -(h * L) - (1.3877787807814457 * 10 ** -16 * K2 * L) / h ** 2 + (
                    h ** 3 * L ** 3) / 10 + 6.938893903907228 * 10 ** -17 * K2 * L ** 3 + (
                                 11 * h ** 5 * L ** 5) / 1500 - (h ** 2 * K2 * L ** 5) / 20 - (
                                 44 * h ** 7 * L ** 7) / 39375 + (h ** 4 * K2 * L ** 7) / 175 + (
                                 h ** 9 * L ** 9) / 17500 - (h ** 6 * K2 * L ** 9) / 3500 + (
                                 8 * h ** 8 * K2 * L ** 11) / 928125 - 1.7878417878417885 * 10 ** -7 * h ** 10 * K2 * L ** 13
        T[1, 2, 2] = -(h ** 3 * L) / 10 + K2 * L + (h ** 5 * L ** 3) / 75 - (h ** 2 * K2 * L ** 3) / 5 - (
                    h ** 7 * L ** 5) / 1875 + (4 * h ** 4 * K2 * L ** 5) / 375 + (2 * h ** 9 * L ** 7) / 196875 - (
                                 2 * h ** 6 * K2 * L ** 7) / 7875 + (2 * h ** 8 * K2 * L ** 9) / 590625 - (
                                 2 * h ** 10 * K2 * L ** 11) / 69609375
        T[1, 2, 3] = K2 * L ** 2 - (2 * h ** 2 * K2 * L ** 4) / 15 + (2 * h ** 4 * K2 * L ** 6) / 375 - (
                    4 * h ** 6 * K2 * L ** 8) / 39375 + (2 * h ** 8 * K2 * L ** 10) / 1771875 - (
                                 4 * h ** 10 * K2 * L ** 12) / 487265625
        T[1, 3, 3] = -(h * L) / 2 + (h ** 3 * L ** 3) / 15 + (K2 * L ** 3) / 3 - (h ** 5 * L ** 5) / 375 - (
                    2 * h ** 2 * K2 * L ** 5) / 75 + (2 * h ** 7 * L ** 7) / 39375 + (
                                 2 * h ** 4 * K2 * L ** 7) / 2625 - (h ** 9 * L ** 9) / 1771875 - (
                                 4 * h ** 6 * K2 * L ** 9) / 354375 + (2 * h ** 8 * K2 * L ** 11) / 19490625 - (
                                 4 * h ** 10 * K2 * L ** 13) / 6334453125
        T[2, 0, 2] = -(h ** 3 * L ** 2) / 5 + K2 * L ** 2 + (h ** 5 * L ** 4) / 30 - (h ** 2 * K2 * L ** 4) / 10 - (
                    91 * h ** 7 * L ** 6) / 45000 + (47 * h ** 4 * K2 * L ** 6) / 9000 + (
                                 41 * h ** 9 * L ** 8) / 630000 - (103 * h ** 6 * K2 * L ** 8) / 630000 + (
                                 1231 * h ** 8 * K2 * L ** 10) / 378000000 - (
                                 16609 * h ** 10 * K2 * L ** 12) / 374220000000
        T[2, 0, 3] = h * L - (7 * h ** 3 * L ** 3) / 30 + (K2 * L ** 3) / 3 + (61 * h ** 5 * L ** 5) / 3000 - (
                    7 * h ** 2 * K2 * L ** 5) / 150 - (547 * h ** 7 * L ** 7) / 630000 + (
                                 3 * h ** 4 * K2 * L ** 7) / 1400 + (703 * h ** 9 * L ** 9) / 32400000 - (
                                 307 * h ** 6 * K2 * L ** 9) / 5670000 + (
                                 11069 * h ** 8 * K2 * L ** 11) / 12474000000 - (
                                 16607 * h ** 10 * K2 * L ** 13) / 1621620000000
        T[2, 1, 2] = -(h ** 3 * L ** 3) / 10 + (K2 * L ** 3) / 3 + (h ** 5 * L ** 5) / 100 - (
                    2 * h ** 2 * K2 * L ** 5) / 75 - (13 * h ** 7 * L ** 7) / 30000 + (
                                 23 * h ** 4 * K2 * L ** 7) / 21000 + (41 * h ** 9 * L ** 9) / 3780000 - (
                                 11 * h ** 6 * K2 * L ** 9) / 405000 + (791 * h ** 8 * K2 * L ** 11) / 1782000000 - (
                                 173 * h ** 10 * K2 * L ** 13) / 33783750000
        T[2, 1, 3] = (h * L ** 2) / 2 - (h ** 3 * L ** 4) / 12 + (K2 * L ** 4) / 6 + (91 * h ** 5 * L ** 6) / 18000 - (
                    11 * h ** 2 * K2 * L ** 6) / 900 - (41 * h ** 7 * L ** 8) / 252000 + (
                                 17 * h ** 4 * K2 * L ** 8) / 42000 + (7381 * h ** 9 * L ** 10) / 2268000000 - (
                                 461 * h ** 6 * K2 * L ** 10) / 56700000 + (
                                 8303 * h ** 8 * K2 * L ** 12) / 74844000000 - (
                                 24911 * h ** 10 * K2 * L ** 14) / 22702680000000
        T[2, 2, 5] = (h ** 2 * L ** 2) / 10 + (1.1102230246251565 * 10 ** -16 * K2 * L ** 2) / h - (
                    11 * h ** 4 * L ** 4) / 300 + (h * K2 * L ** 4) / 12 + (223 * h ** 6 * L ** 6) / 90000 - (
                                 11 * h ** 3 * K2 * L ** 6) / 1800 - (73 * h ** 8 * L ** 8) / 900000 + (
                                 17 * h ** 5 * K2 * L ** 8) / 84000 + 1.6265432098765437 * 10 ** -6 * h ** 10 * L ** 10 - (
                                 461 * h ** 7 * K2 * L ** 10) / 113400000 + 5.546870824648605 * 10 ** -8 * h ** 9 * K2 * L ** 12
        T[2, 3, 5] = (h ** 2 * L ** 3) / 5 - (3 * h ** 4 * L ** 5) / 125 + (h * K2 * L ** 5) / 20 + (
                    113 * h ** 6 * L ** 7) / 105000 - (11 * h ** 3 * K2 * L ** 7) / 4200 - (
                                 16 * h ** 8 * L ** 9) / 590625 + (17 * h ** 5 * K2 * L ** 9) / 252000 + (
                                 41 * h ** 10 * L ** 11) / 92400000 - (461 * h ** 7 * K2 * L ** 11) / 415800000 + (
                                 8303 * h ** 9 * K2 * L ** 13) / 648648000000
        T[3, 0, 2] = -(h ** 3 * L) / 5 + 2 * K2 * L + (7 * h ** 5 * L ** 3) / 150 - (2 * h ** 2 * K2 * L ** 3) / 5 - (
                    61 * h ** 7 * L ** 5) / 15000 + (47 * h ** 4 * K2 * L ** 5) / 1500 + (
                                 547 * h ** 9 * L ** 7) / 3150000 - (103 * h ** 6 * K2 * L ** 7) / 78750 + (
                                 1231 * h ** 8 * K2 * L ** 9) / 37800000 - (
                                 16609 * h ** 10 * K2 * L ** 11) / 31185000000
        T[3, 0, 3] = -(h ** 3 * L ** 2) / 5 + K2 * L ** 2 + (h ** 5 * L ** 4) / 30 - (7 * h ** 2 * K2 * L ** 4) / 30 - (
                    91 * h ** 7 * L ** 6) / 45000 + (3 * h ** 4 * K2 * L ** 6) / 200 + (
                                 41 * h ** 9 * L ** 8) / 630000 - (307 * h ** 6 * K2 * L ** 8) / 630000 + (
                                 11069 * h ** 8 * K2 * L ** 10) / 1134000000 - (
                                 16607 * h ** 10 * K2 * L ** 12) / 124740000000
        T[3, 1, 2] = -(h ** 3 * L ** 2) / 10 + K2 * L ** 2 + (h ** 5 * L ** 4) / 60 - (
                    2 * h ** 2 * K2 * L ** 4) / 15 - (91 * h ** 7 * L ** 6) / 90000 + (
                                 23 * h ** 4 * K2 * L ** 6) / 3000 + (41 * h ** 9 * L ** 8) / 1260000 - (
                                 11 * h ** 6 * K2 * L ** 8) / 45000 + (791 * h ** 8 * K2 * L ** 10) / 162000000 - (
                                 173 * h ** 10 * K2 * L ** 12) / 2598750000
        T[3, 1, 3] = -(h ** 3 * L ** 3) / 10 + (2 * K2 * L ** 3) / 3 + (h ** 5 * L ** 5) / 100 - (
                    11 * h ** 2 * K2 * L ** 5) / 150 - (13 * h ** 7 * L ** 7) / 30000 + (
                                 17 * h ** 4 * K2 * L ** 7) / 5250 + (41 * h ** 9 * L ** 9) / 3780000 - (
                                 461 * h ** 6 * K2 * L ** 9) / 5670000 + (8303 * h ** 8 * K2 * L ** 11) / 6237000000 - (
                                 24911 * h ** 10 * K2 * L ** 13) / 1621620000000
        T[3, 2, 5] = (h ** 2 * L) / 5 - (7 * h ** 4 * L ** 3) / 150 + (h * K2 * L ** 3) / 3 + (
                    73 * h ** 6 * L ** 5) / 15000 - (11 * h ** 3 * K2 * L ** 5) / 300 - (
                                 97 * h ** 8 * L ** 7) / 450000 + (17 * h ** 5 * K2 * L ** 7) / 10500 + (
                                 1229 * h ** 10 * L ** 9) / 226800000 - (461 * h ** 7 * K2 * L ** 9) / 11340000 + (
                                 8303 * h ** 9 * K2 * L ** 11) / 12474000000
        T[3, 3, 5] = (h ** 2 * L ** 2) / 10 - (2.220446049250313 * 10 ** -16 * K2 * L ** 2) / h - (
                    11 * h ** 4 * L ** 4) / 300 + (h * K2 * L ** 4) / 4 + (223 * h ** 6 * L ** 6) / 90000 - (
                                 11 * h ** 3 * K2 * L ** 6) / 600 - (73 * h ** 8 * L ** 8) / 900000 + (
                                 17 * h ** 5 * K2 * L ** 8) / 28000 + 1.6265432098765437 * 10 ** -6 * h ** 10 * L ** 10 - (
                                 461 * h ** 7 * K2 * L ** 10) / 37800000 + 1.6640612473945816 * 10 ** -7 * h ** 9 * K2 * L ** 12
        return T

    if np.abs(kx2 * ky2) < 1e-3:

        T[0, 0, 0] = -(h ** 3 * L ** 2) / 2 - (K2 * L ** 2) / 2 + (h ** 5 * L ** 4) / 6 + (h ** 2 * K2 * L ** 4) / 8 - (
                    h ** 7 * L ** 6) / 45 - (11 * h ** 4 * K2 * L ** 6) / 720 + (h ** 9 * L ** 8) / 630 + (
                                 43 * h ** 6 * K2 * L ** 8) / 40320 - (
                                 19 * h ** 8 * K2 * L ** 10) / 403200 + 1.4258825022713913 * 10 ** -6 * h ** 10 * K2 * L ** 12
        T[0, 0, 1] = h * L - (2 * h ** 3 * L ** 3) / 3 - (K2 * L ** 3) / 3 + (2 * h ** 5 * L ** 5) / 15 + (
                    h ** 2 * K2 * L ** 5) / 12 - (4 * h ** 7 * L ** 7) / 315 - (h ** 4 * K2 * L ** 7) / 120 + (
                                 2 * h ** 9 * L ** 9) / 2835 + (17 * h ** 6 * K2 * L ** 9) / 36288 - (
                                 31 * h ** 8 * K2 * L ** 11) / 1814400 + (h ** 10 * K2 * L ** 13) / 2280960
        T[0, 0, 5] = h ** 2 * L ** 2 - (h ** 4 * L ** 4) / 3 - (h * K2 * L ** 4) / 12 + (2 * h ** 6 * L ** 6) / 45 + (
                    h ** 3 * K2 * L ** 6) / 45 - (h ** 8 * L ** 8) / 315 - (13 * h ** 5 * K2 * L ** 8) / 6720 + (
                                 2 * h ** 10 * L ** 10) / 14175 + (83 * h ** 7 * K2 * L ** 10) / 907200 - (
                                 677 * h ** 9 * K2 * L ** 12) / 239500800
        T[0, 1, 1] = (h * L ** 2) / 4 - (7 * h ** 3 * L ** 4) / 48 - (K2 * L ** 4) / 12 + (
                    31 * h ** 5 * L ** 6) / 1440 + (h ** 2 * K2 * L ** 6) / 72 - (127 * h ** 7 * L ** 8) / 80640 - (
                                 h ** 4 * K2 * L ** 8) / 960 + (73 * h ** 9 * L ** 10) / 1036800 + (
                                 17 * h ** 6 * K2 * L ** 10) / 362880 - (31 * h ** 8 * K2 * L ** 12) / 21772800 + (
                                 h ** 10 * K2 * L ** 14) / 31933440
        T[0, 1, 5] = (h ** 2 * L ** 3) / 2 - (h ** 4 * L ** 5) / 8 - (h * K2 * L ** 5) / 20 + (h ** 6 * L ** 7) / 80 + (
                    h ** 3 * K2 * L ** 7) / 140 - (17 * h ** 8 * L ** 9) / 24192 - (h ** 5 * K2 * L ** 9) / 2240 + (
                                 31 * h ** 10 * L ** 11) / 1209600 + (h ** 7 * K2 * L ** 11) / 59400 - (
                                 151 * h ** 9 * K2 * L ** 13) / 345945600
        T[0, 5, 5] = -(h * L ** 2) / 2  + (h ** 3 * L ** 4) / 6 - (h ** 5 * L ** 6) / 45 - (h ** 2 * K2 * L ** 6) / 120 + (
                                 h ** 7 * L ** 8) / 630 + (h ** 4 * K2 * L ** 8) / 1120 - (h ** 9 * L ** 10) / 14175 - (
                                 h ** 6 * K2 * L ** 10) / 22400 + (
                                 h ** 8 * K2 * L ** 12) / 712800 - 3.117748653462939 * 10 ** -8 * h ** 10 * K2 * L ** 14
        T[0, 2, 2] = (K2 * L ** 2) / 2 - (h ** 2 * K2 * L ** 4) / 24 + (h ** 4 * K2 * L ** 6) / 720 - (
                    h ** 6 * K2 * L ** 8) / 40320 + (h ** 8 * K2 * L ** 10) / 3628800 - (
                                 h ** 10 * K2 * L ** 12) / 479001600
        T[0, 2, 3] = (K2 * L ** 3) / 3 - (h ** 2 * K2 * L ** 5) / 60 + (h ** 4 * K2 * L ** 7) / 2520 - (
                    h ** 6 * K2 * L ** 9) / 181440 + (h ** 8 * K2 * L ** 11) / 19958400 - (
                                 h ** 10 * K2 * L ** 13) / 3113510400
        T[0, 3, 3] = -(h * L ** 2) / 4 + (h ** 3 * L ** 4) / 48 + (K2 * L ** 4) / 12 - (h ** 5 * L ** 6) / 1440 - (
                    h ** 2 * K2 * L ** 6) / 360 + (h ** 7 * L ** 8) / 80640 + (h ** 4 * K2 * L ** 8) / 20160 - (
                                 h ** 9 * L ** 10) / 7257600 - (h ** 6 * K2 * L ** 10) / 1814400 + (
                                 h ** 8 * K2 * L ** 12) / 239500800 - (h ** 10 * K2 * L ** 14) / 43589145600
        T[1, 0, 0] = 7.401486830834378 * 10 ** -17 * h ** 3 * L - K2 * L - 4.934324553889585 * 10 ** -17 * h ** 5 * L ** 3 + (
                    h ** 2 * K2 * L ** 3) / 2 + 9.86864910777917 * 10 ** -18 * h ** 7 * L ** 5 - (
                                   11 * h ** 4 * K2 * L ** 5) / 120 - 9.398713435980162 * 10 ** -19 * h ** 9 * L ** 7 + (
                                   43 * h ** 6 * K2 * L ** 7) / 5040 - (19 * h ** 8 * K2 * L ** 9) / 40320 + (
                                   683 * h ** 10 * K2 * L ** 11) / 39916800
        T[1, 0, 1] = -(K2 * L ** 2) + (5 * h ** 2 * K2 * L ** 4) / 12 - (7 * h ** 4 * K2 * L ** 6) / 120 + (
                    17 * h ** 6 * K2 * L ** 8) / 4032 - (341 * h ** 8 * K2 * L ** 10) / 1814400 + (
                                 13 * h ** 10 * K2 * L ** 12) / 2280960
        T[1, 0, 5] = h ** 2 * L - (h ** 4 * L ** 3) / 6 - (h * K2 * L ** 3) / 3 + (h ** 6 * L ** 5) / 120 + (
                    2 * h ** 3 * K2 * L ** 5) / 15 - (h ** 8 * L ** 7) / 5040 - (13 * h ** 5 * K2 * L ** 7) / 840 + (
                                 h ** 10 * L ** 9) / 362880 + (83 * h ** 7 * K2 * L ** 9) / 90720 - (
                                 677 * h ** 9 * K2 * L ** 11) / 19958400
        T[1, 1, 1] = -(h * L) / 2 + (h ** 3 * L ** 3) / 12 - (K2 * L ** 3) / 3 - (h ** 5 * L ** 5) / 240 + (
                    h ** 2 * K2 * L ** 5) / 12 + (h ** 7 * L ** 7) / 10080 - (h ** 4 * K2 * L ** 7) / 120 - (
                                 h ** 9 * L ** 9) / 725760 + (17 * h ** 6 * K2 * L ** 9) / 36288 - (
                                 31 * h ** 8 * K2 * L ** 11) / 1814400 + (h ** 10 * K2 * L ** 13) / 2280960
        T[1, 1, 5] = -(h * K2 * L ** 4) / 4 + (h ** 3 * K2 * L ** 6) / 20 - (9 * h ** 5 * K2 * L ** 8) / 2240 + (
                    h ** 7 * K2 * L ** 10) / 5400 - (151 * h ** 9 * K2 * L ** 12) / 26611200
        T[1, 5, 5] = -(h * L) +  (h ** 3 * L ** 3) / 6 - (h ** 5 * L ** 5) / 120 - (h ** 2 * K2 * L ** 5) / 20 + \
                     (h ** 7 * L ** 7) / 5040 + (h ** 4 * K2 * L ** 7) / 140 - (h ** 9 * L ** 9) / 362880 - (
                                 h ** 6 * K2 * L ** 9) / 2240 + (h ** 8 * K2 * L ** 11) / 59400 - (
                                 151 * h ** 10 * K2 * L ** 13) / 345945600
        T[1, 2, 2] = K2 * L - (h ** 2 * K2 * L ** 3) / 6 + (h ** 4 * K2 * L ** 5) / 120 - (
                    h ** 6 * K2 * L ** 7) / 5040 + (h ** 8 * K2 * L ** 9) / 362880 - (h ** 10 * K2 * L ** 11) / 39916800
        T[1, 2, 3] = K2 * L ** 2 - (h ** 2 * K2 * L ** 4) / 12 + (h ** 4 * K2 * L ** 6) / 360 - (
                    h ** 6 * K2 * L ** 8) / 20160 + (h ** 8 * K2 * L ** 10) / 1814400 - (
                                 h ** 10 * K2 * L ** 12) / 239500800
        T[1, 3, 3] = -(h * L) / 2 + (h ** 3 * L ** 3) / 12 + (K2 * L ** 3) / 3 - (h ** 5 * L ** 5) / 240 - (
                    h ** 2 * K2 * L ** 5) / 60 + (h ** 7 * L ** 7) / 10080 + (h ** 4 * K2 * L ** 7) / 2520 - (
                                 h ** 9 * L ** 9) / 725760 - (h ** 6 * K2 * L ** 9) / 181440 + (
                                 h ** 8 * K2 * L ** 11) / 19958400 - (h ** 10 * K2 * L ** 13) / 3113510400
        T[2, 0, 2] = K2 * L ** 2 - (h ** 2 * K2 * L ** 4) / 12 + (h ** 4 * K2 * L ** 6) / 360 - (
                    h ** 6 * K2 * L ** 8) / 20160 + (h ** 8 * K2 * L ** 10) / 1814400 - (
                                 h ** 10 * K2 * L ** 12) / 239500800
        T[2, 0, 3] = h * L - (h ** 3 * L ** 3) / 6 + (K2 * L ** 3) / 3 + (h ** 5 * L ** 5) / 120 - (
                    h ** 2 * K2 * L ** 5) / 20 - (h ** 7 * L ** 7) / 5040 + (h ** 4 * K2 * L ** 7) / 504 + (
                                 h ** 9 * L ** 9) / 362880 - (h ** 6 * K2 * L ** 9) / 25920 + (
                                 h ** 8 * K2 * L ** 11) / 2217600 - (h ** 10 * K2 * L ** 13) / 283046400
        T[2, 1, 2] = (K2 * L ** 3) / 3 - (h ** 2 * K2 * L ** 5) / 60 + (h ** 4 * K2 * L ** 7) / 2520 - (
                    h ** 6 * K2 * L ** 9) / 181440 + (h ** 8 * K2 * L ** 11) / 19958400 - (
                                 h ** 10 * K2 * L ** 13) / 3113510400
        T[2, 1, 3] = (h * L ** 2) / 2 - (h ** 3 * L ** 4) / 24 + (K2 * L ** 4) / 6 + (h ** 5 * L ** 6) / 720 - (
                    h ** 2 * K2 * L ** 6) / 90 - (h ** 7 * L ** 8) / 40320 + (h ** 4 * K2 * L ** 8) / 3360 + (
                                 h ** 9 * L ** 10) / 3628800 - (h ** 6 * K2 * L ** 10) / 226800 + (
                                 h ** 8 * K2 * L ** 12) / 23950080 - (h ** 10 * K2 * L ** 14) / 3632428800
        T[2, 2, 5] = (h * K2 * L ** 4) / 12 - (h ** 3 * K2 * L ** 6) / 360 + (h ** 5 * K2 * L ** 8) / 20160 - (
                    h ** 7 * K2 * L ** 10) / 1814400 + (h ** 9 * K2 * L ** 12) / 239500800
        T[2, 3, 5] = (h ** 2 * L ** 3) / 6 - (h ** 4 * L ** 5) / 120 + (h * K2 * L ** 5) / 20 + (
                    h ** 6 * L ** 7) / 5040 - (h ** 3 * K2 * L ** 7) / 504 - (h ** 8 * L ** 9) / 362880 + (
                                 h ** 5 * K2 * L ** 9) / 25920 + (h ** 10 * L ** 11) / 39916800 - (
                                 h ** 7 * K2 * L ** 11) / 2217600 + (h ** 9 * K2 * L ** 13) / 283046400
        T[3, 0, 2] = 2 * K2 * L - (h ** 2 * K2 * L ** 3) / 3 + (h ** 4 * K2 * L ** 5) / 60 - (
                    h ** 6 * K2 * L ** 7) / 2520 + (h ** 8 * K2 * L ** 9) / 181440 - (h ** 10 * K2 * L ** 11) / 19958400
        T[3, 0, 3] = K2 * L ** 2 - (h ** 2 * K2 * L ** 4) / 4 + (h ** 4 * K2 * L ** 6) / 72 - (
                    h ** 6 * K2 * L ** 8) / 2880 + (h ** 8 * K2 * L ** 10) / 201600 - (
                                 h ** 10 * K2 * L ** 12) / 21772800
        T[3, 1, 2] = K2 * L ** 2 - (h ** 2 * K2 * L ** 4) / 12 + (h ** 4 * K2 * L ** 6) / 360 - (
                    h ** 6 * K2 * L ** 8) / 20160 + (h ** 8 * K2 * L ** 10) / 1814400 - (
                                 h ** 10 * K2 * L ** 12) / 239500800
        T[3, 1, 3] = (2 * K2 * L ** 3) / 3 - (h ** 2 * K2 * L ** 5) / 15 + (h ** 4 * K2 * L ** 7) / 420 - (
                    h ** 6 * K2 * L ** 9) / 22680 + (h ** 8 * K2 * L ** 11) / 1995840 - (
                                 h ** 10 * K2 * L ** 13) / 259459200
        T[3, 2, 5] = (h * K2 * L ** 3) / 3 - (h ** 3 * K2 * L ** 5) / 60 + (h ** 5 * K2 * L ** 7) / 2520 - (
                    h ** 7 * K2 * L ** 9) / 181440 + (h ** 9 * K2 * L ** 11) / 19958400
        T[3, 3, 5] = (h * K2 * L ** 4) / 4 - (h ** 3 * K2 * L ** 6) / 72 + (h ** 5 * K2 * L ** 8) / 2880 - (
                    h ** 7 * K2 * L ** 10) / 201600 + (h ** 9 * K2 * L ** 12) / 21772800
        return T

    if h < 1e-6:

        T[0, 0, 0] = -(K2 * L ** 2) / 2 + (K1 * K2 * L ** 4) / 8 - (11 * K1 ** 2 * K2 * L ** 6) / 720 + (
                43 * K1 ** 3 * K2 * L ** 8) / 40320 - (19 * K1 ** 4 * K2 * L ** 10) / 403200 + (
                             683 * K1 ** 5 * K2 * L ** 12) / 479001600 - (
                             2731 * K1 ** 6 * K2 * L ** 14) / 87178291200
        T[0, 0, 1] = -(K2 * L ** 3) / 3 + (K1 * K2 * L ** 5) / 12 - (K1 ** 2 * K2 * L ** 7) / 120 + (
                17 * K1 ** 3 * K2 * L ** 9) / 36288 - (31 * K1 ** 4 * K2 * L ** 11) / 1814400 + (
                             K1 ** 5 * K2 * L ** 13) / 2280960 - (5461 * K1 ** 6 * K2 * L ** 15) / 653837184000
        T[0, 0, 5] = (K1 * L ** 2) / 2 - (K1 ** 2 * L ** 4) / 12 + (K1 ** 3 * L ** 6) / 240 - (
                K1 ** 4 * L ** 8) / 10080 + (K1 ** 5 * L ** 10) / 725760 - (K1 ** 6 * L ** 12) / 79833600
        T[0, 1, 1] = -(K2 * L ** 4) / 12 + (K1 * K2 * L ** 6) / 72 - (K1 ** 2 * K2 * L ** 8) / 960 + (
                17 * K1 ** 3 * K2 * L ** 10) / 362880 - (31 * K1 ** 4 * K2 * L ** 12) / 21772800 + (
                             K1 ** 5 * K2 * L ** 14) / 31933440 - (5461 * K1 ** 6 * K2 * L ** 16) / 10461394944000
        T[0, 1, 5] = (K1 * L ** 3) / 6 - (K1 ** 2 * L ** 5) / 60 + (K1 ** 3 * L ** 7) / 1680 - (
                K1 ** 4 * L ** 9) / 90720 + (K1 ** 5 * L ** 11) / 7983360 - (K1 ** 6 * L ** 13) / 1037836800
        T[0, 5, 5] = 0
        T[0, 2, 2] = (K2 * L ** 2) / 2 + (K1 * K2 * L ** 4) / 24 + (7 * K1 ** 2 * K2 * L ** 6) / 720 + (
                5 * K1 ** 3 * K2 * L ** 8) / 8064 + (103 * K1 ** 4 * K2 * L ** 10) / 3628800 + (
                             409 * K1 ** 5 * K2 * L ** 12) / 479001600 + (149 * K1 ** 6 * K2 * L ** 14) / 7925299200
        T[0, 2, 3] = (K2 * L ** 3) / 3 + (K1 * K2 * L ** 5) / 20 + (13 * K1 ** 2 * K2 * L ** 7) / 2520 + (
                17 * K1 ** 3 * K2 * L ** 9) / 60480 + (41 * K1 ** 4 * K2 * L ** 11) / 3991680 + (
                             K1 ** 5 * K2 * L ** 13) / 3801600 + (3277 * K1 ** 6 * K2 * L ** 15) / 653837184000
        T[0, 3, 3] = (K2 * L ** 4) / 12 + (K1 * K2 * L ** 6) / 120 + (13 * K1 ** 2 * K2 * L ** 8) / 20160 + (
                17 * K1 ** 3 * K2 * L ** 10) / 604800 + (41 * K1 ** 4 * K2 * L ** 12) / 47900160 + (
                             K1 ** 5 * K2 * L ** 14) / 53222400 + (3277 * K1 ** 6 * K2 * L ** 16) / 10461394944000
        T[1, 0, 0] = -(K2 * L) + (K1 * K2 * L ** 3) / 2 - (11 * K1 ** 2 * K2 * L ** 5) / 120 + (
                43 * K1 ** 3 * K2 * L ** 7) / 5040 - (19 * K1 ** 4 * K2 * L ** 9) / 40320 + (
                             683 * K1 ** 5 * K2 * L ** 11) / 39916800 - (2731 * K1 ** 6 * K2 * L ** 13) / 6227020800
        T[1, 0, 1] = -(K2 * L ** 2) + (5 * K1 * K2 * L ** 4) / 12 - (7 * K1 ** 2 * K2 * L ** 6) / 120 + (
                17 * K1 ** 3 * K2 * L ** 8) / 4032 - (341 * K1 ** 4 * K2 * L ** 10) / 1814400 + (
                             13 * K1 ** 5 * K2 * L ** 12) / 2280960 - (5461 * K1 ** 6 * K2 * L ** 14) / 43589145600
        T[1, 0, 5] = K1 * L - (K1 ** 2 * L ** 3) / 3 + (K1 ** 3 * L ** 5) / 40 - (K1 ** 4 * L ** 7) / 1260 + (
                K1 ** 5 * L ** 9) / 72576 - (K1 ** 6 * L ** 11) / 6652800
        T[1, 1, 1] = -(K2 * L ** 3) / 3 + (K1 * K2 * L ** 5) / 12 - (K1 ** 2 * K2 * L ** 7) / 120 + (
                17 * K1 ** 3 * K2 * L ** 9) / 36288 - (31 * K1 ** 4 * K2 * L ** 11) / 1814400 + (
                             K1 ** 5 * K2 * L ** 13) / 2280960 - (5461 * K1 ** 6 * K2 * L ** 15) / 653837184000
        T[1, 1, 5] = (K1 * L ** 2) / 2 - (K1 ** 2 * L ** 4) / 12 + (K1 ** 3 * L ** 6) / 240 - (
                K1 ** 4 * L ** 8) / 10080 + (K1 ** 5 * L ** 10) / 725760 - (K1 ** 6 * L ** 12) / 79833600
        T[1, 5, 5] = 0
        T[1, 2, 2] = K2 * L + (K1 * K2 * L ** 3) / 6 + (7 * K1 ** 2 * K2 * L ** 5) / 120 + (
                5 * K1 ** 3 * K2 * L ** 7) / 1008 + (103 * K1 ** 4 * K2 * L ** 9) / 362880 + (
                             409 * K1 ** 5 * K2 * L ** 11) / 39916800 + (149 * K1 ** 6 * K2 * L ** 13) / 566092800
        T[1, 2, 3] = K2 * L ** 2 + (K1 * K2 * L ** 4) / 4 + (13 * K1 ** 2 * K2 * L ** 6) / 360 + (
                17 * K1 ** 3 * K2 * L ** 8) / 6720 + (41 * K1 ** 4 * K2 * L ** 10) / 362880 + (
                             13 * K1 ** 5 * K2 * L ** 12) / 3801600 + (3277 * K1 ** 6 * K2 * L ** 14) / 43589145600
        T[1, 3, 3] = (K2 * L ** 3) / 3 + (K1 * K2 * L ** 5) / 20 + (13 * K1 ** 2 * K2 * L ** 7) / 2520 + (
                17 * K1 ** 3 * K2 * L ** 9) / 60480 + (41 * K1 ** 4 * K2 * L ** 11) / 3991680 + (
                             K1 ** 5 * K2 * L ** 13) / 3801600 + (3277 * K1 ** 6 * K2 * L ** 15) / 653837184000
        T[2, 0, 2] = K2 * L ** 2 + (K1 * K2 * L ** 4) / 12 - (K1 ** 2 * K2 * L ** 6) / 120 - (
                K1 ** 3 * K2 * L ** 8) / 6720 + (13 * K1 ** 4 * K2 * L ** 10) / 1814400 + (
                             13 * K1 ** 5 * K2 * L ** 12) / 239500800 - (17 * K1 ** 6 * K2 * L ** 14) / 14529715200
        T[2, 0, 3] = (K2 * L ** 3) / 3 - (K1 * K2 * L ** 5) / 60 - (K1 ** 2 * K2 * L ** 7) / 504 + (
                K1 ** 3 * K2 * L ** 9) / 60480 + (19 * K1 ** 4 * K2 * L ** 11) / 19958400 - (
                             K1 ** 5 * K2 * L ** 13) / 239500800 - (K1 ** 6 * K2 * L ** 15) / 8491392000
        T[2, 1, 2] = (K2 * L ** 3) / 3 + (K1 * K2 * L ** 5) / 20 - (K1 ** 2 * K2 * L ** 7) / 2520 - (
                K1 ** 3 * K2 * L ** 9) / 20160 + (K1 ** 4 * K2 * L ** 11) / 2851200 + (
                             K1 ** 5 * K2 * L ** 13) / 79833600 - (K1 ** 6 * K2 * L ** 15) / 26153487360
        T[2, 1, 3] = (K2 * L ** 4) / 6 + (K1 * K2 * L ** 6) / 180 - (K1 ** 2 * K2 * L ** 8) / 3360 - (
                K1 ** 3 * K2 * L ** 10) / 302400 + (13 * K1 ** 4 * K2 * L ** 12) / 119750400 + (
                             K1 ** 5 * K2 * L ** 14) / 1676505600 - (17 * K1 ** 6 * K2 * L ** 16) / 1743565824000
        T[2, 2, 5] = -(K1 * L ** 2) / 2 - (K1 ** 2 * L ** 4) / 12 - (K1 ** 3 * L ** 6) / 240 - (
                K1 ** 4 * L ** 8) / 10080 - (K1 ** 5 * L ** 10) / 725760 - (K1 ** 6 * L ** 12) / 79833600
        T[2, 3, 5] = -(K1 * L ** 3) / 6 - (K1 ** 2 * L ** 5) / 60 - (K1 ** 3 * L ** 7) / 1680 - (
                K1 ** 4 * L ** 9) / 90720 - (K1 ** 5 * L ** 11) / 7983360 - (K1 ** 6 * L ** 13) / 1037836800
        T[3, 0, 2] = 2 * K2 * L + (K1 * K2 * L ** 3) / 3 - (K1 ** 2 * K2 * L ** 5) / 20 - (
                K1 ** 3 * K2 * L ** 7) / 840 + (13 * K1 ** 4 * K2 * L ** 9) / 181440 + (
                             13 * K1 ** 5 * K2 * L ** 11) / 19958400 - (17 * K1 ** 6 * K2 * L ** 13) / 1037836800
        T[3, 0, 3] = K2 * L ** 2 - (K1 * K2 * L ** 4) / 12 - (K1 ** 2 * K2 * L ** 6) / 72 + (
                K1 ** 3 * K2 * L ** 8) / 6720 + (19 * K1 ** 4 * K2 * L ** 10) / 1814400 - (
                             13 * K1 ** 5 * K2 * L ** 12) / 239500800 - (K1 ** 6 * K2 * L ** 14) / 566092800
        T[3, 1, 2] = K2 * L ** 2 + (K1 * K2 * L ** 4) / 4 - (
                K1 ** 2 * K2 * L ** 6) / 360 - (K1 ** 3 * K2 * L ** 8) / 2240 + (
                             K1 ** 4 * K2 * L ** 10) / 259200 + (13 * K1 ** 5 * K2 * L ** 12) / 79833600 - (
                             K1 ** 6 * K2 * L ** 14) / 1743565824
        T[3, 1, 3] = (2 * K2 * L ** 3) / 3 + (K1 * K2 * L ** 5) / 30 - (K1 ** 2 * K2 * L ** 7) / 420 - (
                K1 ** 3 * K2 * L ** 9) / 30240 + (13 * K1 ** 4 * K2 * L ** 11) / 9979200 + (
                             K1 ** 5 * K2 * L ** 13) / 119750400 - 1.5600214012912426 * 10 ** -10 * K1 ** 6 * K2 * L ** 15
        T[3, 2, 5] = -(K1 * L) - (K1 ** 2 * L ** 3) / 3 - (K1 ** 3 * L ** 5) / 40 - (K1 ** 4 * L ** 7) / 1260 - (
                K1 ** 5 * L ** 9) / 72576 - (K1 ** 6 * L ** 11) / 6652800
        T[3, 3, 5] = -(K1 * L ** 2) / 2 - (K1 ** 2 * L ** 4) / 12 - (K1 ** 3 * L ** 6) / 240 - (
                K1 ** 4 * L ** 8) / 10080 - (K1 ** 5 * L ** 10) / 725760 - (K1 ** 6 * L ** 12) / 79833600
        return T
