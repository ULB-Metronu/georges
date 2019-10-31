from numba import njit
import numpy as np
from numpy import cos, sin, cosh, sinh, sqrt


@njit
def compute_transport_quadrupole_matrix(element_parameters: list) -> np.ndarray:
    L: float = element_parameters[0]
    K1: float = element_parameters[1]
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
    R[0, 0] = cx
    R[0, 1] = sx
    R[1, 0] = -(K1 * sx)
    R[1, 1] = cx
    R[2, 2] = cy
    R[2, 3] = sy
    R[3, 2] = K1 * sy
    R[3, 3] = cy
    R[4, 4] = 1
    R[5, 5] = 1

    return R


@njit
def compute_transport_quadrupole_tensor(element_parameters: list) -> np.ndarray:
    L: float = element_parameters[0]
    K1: float = element_parameters[1]
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
    T[0, 0, 5] = (K1 * L * sx) / 2
    T[0, 1, 5] = -(cx * L) / 2 - sx / 2  # -(cx * L) / 2 + sx / 2
    T[1, 0, 5] = (cx * K1 * L) / 2 - (K1 * sx) / 2  # (cx * K1 * L) / 2 + (K1 * sx)
    T[1, 1, 5] = (K1 * L * sx) / 2
    T[2, 2, 5] = -(K1 * L * sy) / 2
    T[2, 3, 5] = (-(cy * L) - sy) / 2  # (-(cy * L) + sy) / 2
    T[3, 2, 5] = -(cy * K1 * L) / 2 + (K1 * sy) / 2  # -(cx * K1 * L) / 2 - (K1 * sy) / 2
    T[3, 3, 5] = -(K1 * L * sy) / 2

    return T
