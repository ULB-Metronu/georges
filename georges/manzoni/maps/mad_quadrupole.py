from numba import njit
import numpy as np
from numpy import cos, sin, cosh, sinh, sqrt


@njit
def compute_mad_quadrupole_matrix(L: float, K1: float, beta: float) -> np.ndarray:
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


@njit
def compute_mad_quadrupole_tensor(L: float, K1: float, beta: float) -> np.ndarray:
    gamma = 1 / sqrt(1 - beta ** 2)
    kx2 = K1
    ky2 = -K1
    T = np.zeros((6, 6, 6))

    # Setup basic notations - Horizontal plane
    if K1 == 0:
        cx = 1
        sx = L
    elif K1 > 0:
        cx = cos(sqrt(K1) * L)
        sx = sin(sqrt(K1) * L) / sqrt(K1)
    else:
        cx = cosh(sqrt(-K1) * L)
        sx = sinh(sqrt(-K1) * L) / sqrt(-K1)

    # Setup basic notations - Vertical plane
    if K1 == 0:
        cy = 1
        sy = L
    elif K1 > 0:
        cy = cosh(sqrt(K1) * L)
        sy = sinh(sqrt(K1) * L) / sqrt(K1)
    else:
        cy = cos(sqrt(-K1) * L)
        sy = sin(sqrt(-K1) * L) / sqrt(-K1)

    # Definition of the tensor elements
    T[0, 0, 5] = (1 / (2 * beta)) * K1 * L * sx
    T[0, 1, 5] = -(1 / (2 * beta)) * (sx + L * cx)
    T[1, 0, 5] = -(K1 / (2 * beta)) * (sx - L * cx)
    T[1, 1, 5] = (K1 / (2 * beta)) * L * sx
    T[2, 2, 5] = (1 / (2 * beta)) * (-K1 * L * sy)
    T[2, 3, 5] = -(1 / (2 * beta)) * (sy + L * cy)
    T[3, 2, 5] = (K1 / (2 * beta)) * (sy - L * cy)
    T[3, 3, 5] = -(K1 / (2 * beta)) * L * sy

    return T
