from numba import njit
import numpy as np
from numpy import cos, sin, cosh, sinh, sqrt


@njit
def compute_mad_quadrupole_matrix(L: float, K1: float, beta: float) -> np.ndarray:
    gamma = 1 / sqrt(1 - beta ** 2)
    kx2 = K1
    ky2 = -K1
    R = np.zeros((6, 6))

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

    # Definition of the matrix elements
    R[0, 0] = cos(sqrt(K1) * L)
    R[0, 1] = sin(sqrt(K1) * L) / sqrt(K1)
    R[1, 0] = -(sqrt(K1) * sin(sqrt(K1) * L))
    R[1, 1] = cos(sqrt(K1) * L)
    R[2, 2] = cos(sqrt(-K1) * L)
    R[2, 3] = sin(sqrt(-K1) * L) / sqrt(-K1)
    R[3, 2] = (K1 * sin(sqrt(-K1) * L)) / sqrt(-K1)
    R[3, 3] = cos(sqrt(-K1) * L)
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
    T[0, 0, 0] = 1

    return T
