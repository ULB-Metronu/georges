from numba import njit
import numpy as np
from numpy import cos, sin, cosh, sinh, sqrt


@njit
def compute_transport_quadrupole_matrix(L: float, K1: float, *args) -> np.ndarray:
    kq2 = K1
    R = np.zeros((6, 6))

    def cq(kq2, L):
        if kq2 > 0:
            cq = cos(sqrt(kq2) * L)
        if kq2 < 0:
            cq = cosh(sqrt(-kq2) * L)
        if kq2 == 0:
            cq = 1
        return cq

    def sq(kq2, L):
        if kq2 > 0:
            sq = sin(sqrt(kq2) * L) / sqrt(kq2)
        if kq2 < 0:
            sq = sinh(sqrt(-kq2) * L) / sqrt(-kq2)
        if kq2 == 0:
            sq = L
        return sq

    R[0, 0] = cq(kq2, L)
    R[0, 1] = sq(kq2, L)
    R[1, 0] = -kq2 * sq(kq2, L)
    R[1, 1] = cq(kq2, L)
    R[2, 2] = cq(-kq2, L)
    R[2, 3] = sq(-kq2, L)
    R[3, 2] = kq2 * sq(-kq2, L)
    R[3, 3] = cq(-kq2, L)
    R[4, 4] = 1
    R[5, 5] = 1

    return R


#@njit
def compute_transport_quadrupole_tensor(L: float, K1: float, *args) -> np.ndarray:
    kq2 = K1
    T = np.zeros((6, 6, 6))

    def cq(kq2, L):
        if kq2 > 0:
            cq = cos(sqrt(kq2) * L)
        if kq2 < 0:
            cq = cosh(sqrt(-kq2) * L)
        if kq2 == 0:
            cq = 1
        return cq

    def sq(kq2, L):
        if kq2 > 0:
            sq = sin(sqrt(kq2) * L) / sqrt(kq2)
        if kq2 < 0:
            sq = sinh(sqrt(-kq2) * L) / sqrt(-kq2)
        if kq2 == 0:
            sq = L
        return sq

    T[0, 0, 5] = (1 / 2) * kq2 * L * sq(kq2, L)
    T[0, 1, 5] = (1 / 2) * sq(kq2, L) - (L / 2) * cq(kq2, L)
    T[1, 0, 5] = (kq2 * L / 2) * cq(kq2, L) - kq2 * sq(kq2, L) / 2
    T[1, 1, 5] = kq2 * L * sq(kq2, L) / 2
    T[2, 2, 5] = -kq2 * L * sq(-kq2, L) / 2
    T[2, 3, 5] = sq(-kq2, L) / 2 - L * cq(-kq2, L) / 2
    T[3, 2, 5] = -kq2 * L / 2 * cq(-kq2, L) - kq2 * sq(-kq2, L) / 2
    T[3, 3, 5] = -kq2 * L * sq(-kq2, L) / 2

    return T
