from numba import njit
import numpy as np
from numpy import cos, sin, cosh, sinh, sqrt


@njit
def compute_transport_combined_dipole_matrix(element_parameters: list) -> np.ndarray:
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    K1: float = element_parameters[2]
    h = alpha / L
    kx2 = h ** 2 + K1
    ky2 = -K1
    R = np.zeros((6, 6))

    # Setup basic notations - Horizontal plane
    if h ** 2 + K1 == 0:
        cx = 1
        sx = L
        dx = (h * L ** 2) / 2
    elif h ** 2 + K1 > 0:
        cx = cos(sqrt(h ** 2 + K1) * L)
        sx = sin(sqrt(h ** 2 + K1) * L) / sqrt(h ** 2 + K1)
        dx = (h * (1 - cos(sqrt(h ** 2 + K1) * L))) / (h ** 2 + K1)
    else:
        cx = cosh(sqrt(-h ** 2 - K1) * L)
        sx = sinh(sqrt(-h ** 2 - K1) * L) / sqrt(-h ** 2 - K1)
        dx = (h * (1 - cosh(sqrt(-h ** 2 - K1) * L))) / (h ** 2 + K1)

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
    R[0, 5] = dx
    R[1, 0] = -(kx2 * sx)
    R[1, 1] = cx
    R[1, 5] = h * sx
    R[2, 2] = cy
    R[2, 3] = sy
    R[3, 2] = -(ky2 * sy)
    R[3, 3] = cy
    R[4, 4] = 1
    R[5, 5] = 1

    return R



@njit(cache=True)
def compute_transport_combined_dipole_tensor(L: float, alpha: float, K1: float, K2: float) -> np.ndarray:
    h = alpha / L
    kx2 = h ** 2 + K1
    ky2 = -K1
    T = np.zeros((6, 6, 6))

    # Setup basic notations - Horizontal plane
    if h ** 2 + K1 == 0:
        cx = 1
        sx = L
        dx = (h * L ** 2) / 2
    elif h ** 2 + K1 > 0:
        cx = cos(sqrt(h ** 2 + K1) * L)
        sx = sin(sqrt(h ** 2 + K1) * L) / sqrt(h ** 2 + K1)
        dx = (h * (1 - cos(sqrt(h ** 2 + K1) * L))) / (h ** 2 + K1)
    else:
        cx = cosh(sqrt(-h ** 2 - K1) * L)
        sx = sinh(sqrt(-h ** 2 - K1) * L) / sqrt(-h ** 2 - K1)
        dx = (h * (1 - cosh(sqrt(-h ** 2 - K1) * L))) / (h ** 2 + K1)

    # Setup basic notations - Vertical plane
    if K1 == 0:
        cy = 1
        sy = L
        I34 = L ** 3 / 6
    elif K1 > 0:
        cy = cosh(sqrt(K1) * L)
        sy = sinh(sqrt(K1) * L) / sqrt(K1)
        I34 = -(-(L * cosh(sqrt(K1) * L)) + sinh(sqrt(K1) * L) / sqrt(K1)) / (2 * K1)
    else:
        cy = cos(sqrt(-K1) * L)
        sy = sin(sqrt(-K1) * L) / sqrt(-K1)
        I34 = -(-(L * cos(sqrt(-K1) * L)) + sin(sqrt(-K1) * L) / sqrt(-K1)) / (2 * K1)

    if h != 0:
        T[0, 0, 0] = (2 * h - 2 * cx * h - (2 * dx * (h ** 3 + 2 * h * K1 + K2)) / h + (
                    (-1 + cx) * h ** 2 * sx ** 2) / dx - 2 * (h ** 3 + 2 * h * K1 + K2) * sx ** 2) / 6
        T[0, 0, 1] = (((2 + cx) * h ** 2 - 2 * dx * h ** 3 - 4 * dx * h * K1 - 2 * dx * K2) * sx) / (3 * h)
        T[0, 0, 5] = (-4 * dx ** 2 * (h ** 3 + 2 * h * K1 + K2) + 2 * dx * h * (
                    -2 * (-1 + cx) * h + 3 * (h ** 3 + 2 * h * K1 + K2) * L * sx - 2 * (
                        h ** 3 + 2 * h * K1 + K2) * sx ** 2) + (-1 + cx) * h * sx * (
                                  3 * K1 * L + 2 * h ** 2 * (3 * L + sx))) / (6 * (-1 + cx) * h)
        T[0, 1, 1] = (4 * dx ** 2 * (h ** 3 + 2 * h * K1 + K2) + (-1 + cx) * h ** 3 * sx ** 2 + dx * h * (
                    (-1 + cx) * h - 2 * (h ** 3 + 2 * h * K1 + K2) * sx ** 2)) / (6 * (-1 + cx) * h ** 2)
        T[0, 1, 5] = (dx * h * ((3 * (2 * h ** 2 + K1) * (cx * L - sx)) / ((-1 + cx) * h ** 2) + 2 * sx + (
                    2 * dx * (h ** 3 + 2 * h * K1 + K2) * (3 * cx * L - sx - 2 * cx * sx)) / (
                                            (-1 + cx) ** 2 * h ** 2))) / 6
        T[0, 5, 5] = -(dx * (6 * (-1 + cx) * dx * K1 + 8 * dx ** 2 * K2 + 2 * dx * h ** 4 * sx * (-3 * L + sx) + h * (
                    6 * (-1 + cx) ** 2 + 16 * dx ** 2 * K1 + 3 * (
                        K1 - cx * K1 - 2 * dx * K2) * L * sx + 2 * dx * K2 * sx ** 2) + 2 * dx * h ** 2 * (
                                         -7 + 7 * cx + 2 * K1 * sx * (-3 * L + sx)) + h ** 3 * (
                                         8 * dx ** 2 - (-1 + cx) * sx * (6 * L + sx)))) / (6 * (-1 + cx) ** 2 * h)
        T[0, 2, 2] = (dx * K1) / 2 + (dx * K2) / h + (dx * K1 * K2 * ((-2 * dx) / h + sy ** 2)) / (
                    h - cx * h + 4 * dx * K1)
        T[0, 2, 3] = (-2 * dx * K2 * (sx - cy * sy)) / (h - cx * h + 4 * dx * K1)
        T[0, 3, 3] = -dx / 2 + (dx * K2 * ((-2 * dx) / h + sy ** 2)) / (h - cx * h + 4 * dx * K1)
        T[1, 0, 0] = -((3 * (-1 + cx) * cx * h ** 2 + dx * (
                    (1 + 2 * cx) * h ** 3 + h * (-1 + cx + 2 * K1 + 4 * cx * K1) + K2 + 2 * cx * K2)) * sx) / (3 * dx)
        T[1, 0, 1] = (-((-1 + cx) * (1 + 3 * cx) * h) + (2 * dx * (h ** 3 + 2 * h * K1 + K2)) / h - (
                    (-1 + cx) * h ** 2 * sx ** 2) / dx - 4 * (h ** 3 + 2 * h * K1 + K2) * sx ** 2) / 3
        T[1, 0, 5] = (h ** 2 * (-2 * (-1 + cx) * sx - 6 * cx * sx + 3 * (2 + K1 / h ** 2) * (cx * L + sx) + (
                    2 * dx * (h ** 3 + 2 * h * K1 + K2) * (3 * cx * L + sx - 4 * cx * sx)) / ((-1 + cx) * h ** 2))) / 6
        T[1, 1, 1] = -(((-1 + 4 * cx) * h ** 2 + 4 * dx * h ** 3 + 8 * dx * h * K1 + 4 * dx * K2) * sx) / (6 * h)
        T[1, 1, 5] = (4 * dx ** 2 * (h ** 3 + 2 * h * K1 + K2) + (-1 + cx) * h * (
                    3 * K1 * L + h ** 2 * (6 * L - 2 * sx)) * sx + 2 * dx * h * (
                                  -((-1 + cx) * (1 + 3 * cx) * h) + 3 * (h ** 3 + 2 * h * K1 + K2) * L * sx - 4 * (
                                      h ** 3 + 2 * h * K1 + K2) * sx ** 2)) / (6 * (-1 + cx) * h)
        T[1, 5, 5] = (h * (
                    -6 * sx - 4 * dx * h * sx + (3 * dx * h * (2 + K1 / h ** 2) * (-(cx * L) + sx)) / (1 - cx) + (
                        2 * dx ** 2 * (h ** 3 + 2 * h * K1 + K2) * (3 * cx * L - sx - 2 * cx * sx)) / (
                                (-1 + cx) ** 2 * h))) / 6
        T[1, 2, 2] = (h * K1 * sx) / 2 + K2 * (sx - (2 * dx * K1 * (sx - cy * sy)) / (h - cx * h + 4 * dx * K1))
        T[1, 2, 3] = (-2 * dx * K2 * (-1 + cx - 2 * K1 * sy ** 2)) / (h - cx * h + 4 * dx * K1)
        T[1, 3, 3] = -(h * sx) / 2 - (2 * dx * K2 * (sx - cy * sy)) / (h - cx * h + 4 * dx * K1)
        T[2, 0, 2] = (2 * (-1 + cx) * cy * dx * K2 - K1 * (
                    -((-1 + cx) * h ** 2) + 4 * dx * h * K1 + 4 * dx * K2) * sx * sy) / ((-1 + cx) * h - 4 * dx * K1)
        T[2, 0, 3] = (cy * ((-1 + cx) * h ** 2 - 4 * dx * h * K1 - 4 * dx * K2) * sx + 2 * (1 + cx) * dx * K2 * sy) / (
                    (-1 + cx) * h - 4 * dx * K1)
        T[2, 1, 2] = (2 * cy * dx * h * K2 * sx - dx * (
                    h * K1 * (h - cx * h + 4 * dx * K1) + 2 * (h + 2 * dx * K1) * K2) * sy) / (
                                 h * ((-1 + cx) * h - 4 * dx * K1))
        T[2, 1, 3] = (dx * (cy * h * ((-1 + cx) * h - 4 * dx * K1) - 4 * cy * dx * K2 + 2 * h * K2 * sx * sy)) / (
                    h * ((-1 + cx) * h - 4 * dx * K1))
        T[2, 2, 5] = -(K1 * L * sy) / 2 + (dx * (
                    2 * (-1 + cx) * cy * dx * K2 - ((-1 + cx) * h - 4 * dx * K1) * (h * K1 + K2) * L * sy - K1 * (
                        -((-1 + cx) * h ** 2) + 4 * dx * h * K1 + 4 * dx * K2) * sx * sy)) / (
                                 (-1 + cx) * ((-1 + cx) * h - 4 * dx * K1))
        T[2, 3, 5] = (-((-1 + cx) * dx * h ** 2 * (2 * I34 * K1 - cy * sx + sy)) + h * (
                    I34 * K1 * (-(-1 + cx) ** 2 + 8 * dx ** 2 * K1) - 2 * (
                        -1 + cx) * dx * I34 * K2 + 4 * dx ** 2 * K1 * (-(cy * sx) + sy)) + 2 * dx * (
                                  2 * I34 * K1 * ((-1 + cx) * K1 + 2 * dx * K2) + dx * K2 * (
                                      -2 * cy * sx + sy + cx * sy))) / ((-1 + cx) * ((-1 + cx) * h - 4 * dx * K1))
        T[3, 0, 2] = (cy * (
                    h * K1 * ((-1 + cx) * h - 4 * dx * K1) + 2 * ((-1 + cx) * h - 2 * dx * K1) * K2) * sx - 2 * (
                                  1 + cx) * dx * K1 * K2 * sy) / ((-1 + cx) * h - 4 * dx * K1)
        T[3, 0, 3] = (-2 * (-1 + cx) * cy * dx * K2 + (
                    h * K1 * ((-1 + cx) * h - 4 * dx * K1) + 2 * ((-1 + cx) * h - 2 * dx * K1) * K2) * sx * sy) / (
                                 (-1 + cx) * h - 4 * dx * K1)
        T[3, 1, 2] = (cy * dx * (h * K1 * ((-1 + cx) * h - 4 * dx * K1) - 2 * (
                    h - cx * h + 2 * dx * K1) * K2) - 2 * dx * h * K1 * K2 * sx * sy) / (
                                 h * ((-1 + cx) * h - 4 * dx * K1))
        T[3, 1, 3] = (-2 * cy * dx * h * K2 * sx + dx * (
                    h * K1 * ((-1 + cx) * h - 4 * dx * K1) + 2 * (cx * h - 2 * dx * K1) * K2) * sy) / (
                                 h * ((-1 + cx) * h - 4 * dx * K1))
        T[3, 2, 5] = -(cy * ((-1 + cx) * h - 4 * dx * K1) * (
                    (-1 + cx + 2 * dx * h) * K1 + 2 * dx * K2) * L + 2 * cy * dx * (
                                   h ** 2 * (K1 - cx * K1) + 4 * dx * K1 * K2 + 2 * h * (
                                       2 * dx * K1 ** 2 + K2 - cx * K2)) * sx + (-1 + cx) * (
                                   (-1 + cx) * h * K1 + 2 * dx * h * K2 + 4 * dx * K1 * (-K1 + dx * K2)) * sy) / (
                                 2 * (-1 + cx) * ((-1 + cx) * h - 4 * dx * K1))
        T[3, 3, 5] = -(4 * (-1 + cx) * cy * dx ** 2 * K2 + ((-1 + cx) * h - 4 * dx * K1) * (
                    (-1 + cx + 2 * dx * h) * K1 + 2 * dx * K2) * L * sy + 2 * dx * (
                                   h ** 2 * (K1 - cx * K1) + 4 * dx * K1 * K2 + 2 * h * (
                                       2 * dx * K1 ** 2 + K2 - cx * K2)) * sx * sy) / (
                                 2 * (-1 + cx) * ((-1 + cx) * h - 4 * dx * K1))

    if h == 0 and K1 != 0:
        T[0, 0, 0] = (K2 * (-1 + cx - kx2 * sx ** 2)) / (3 * kx2)
        T[0, 0, 1] = (2 * (-1 + cx) * K2 * sx) / (3 * kx2)
        T[0, 0, 5] = (K1 * L * sx) / 2
        T[0, 1, 1] = (K2 * (-2 + 2 * cx + kx2 * sx ** 2)) / (3 * kx2 ** 2)
        T[0, 1, 5] = (K1 * (-(cx * L) + sx)) / (2 * kx2)
        T[0, 5, 5] = 0
        T[0, 2, 2] = (K2 * (2 * (-1 + cx) * ky2 - kx2 * (-1 + cx + ky2 * sy ** 2))) / (kx2 ** 2 - 4 * kx2 * ky2)
        T[0, 2, 3] = (-2 * K2 * (sx - cy * sy)) / (kx2 - 4 * ky2)
        T[0, 3, 3] = (K2 * (-2 + 2 * cx + kx2 * sy ** 2)) / (kx2 ** 2 - 4 * kx2 * ky2)
        T[1, 0, 0] = -((1 + 2 * cx) * K2 * sx) / 3
        T[1, 0, 1] = (-2 * K2 * (-1 + cx + 2 * kx2 * sx ** 2)) / (3 * kx2)
        T[1, 0, 5] = (K1 * (cx * L + sx)) / 2
        T[1, 1, 1] = (2 * (-1 + cx) * K2 * sx) / (3 * kx2)
        T[1, 1, 5] = (K1 * L * sx) / 2
        T[1, 5, 5] = 0
        T[1, 2, 2] = K2 * (sx - (2 * ky2 * (-sx + cy * sy)) / (kx2 - 4 * ky2))
        T[1, 2, 3] = (-2 * K2 * (-1 + cx + 2 * ky2 * sy ** 2)) / (kx2 - 4 * ky2)
        T[1, 3, 3] = (-2 * K2 * (sx - cy * sy)) / (kx2 - 4 * ky2)
        T[2, 0, 2] = (2 * K2 * (cy - cx * cy - 2 * ky2 * sx * sy)) / (kx2 - 4 * ky2)
        T[2, 0, 3] = (-2 * K2 * (-2 * cy * sx + sy + cx * sy)) / (kx2 - 4 * ky2)
        T[2, 1, 2] = (2 * K2 * (-(cy * kx2 * sx) + (kx2 + 2 * (-1 + cx) * ky2) * sy)) / (kx2 ** 2 - 4 * kx2 * ky2)
        T[2, 1, 3] = (-2 * K2 * (2 * (-1 + cx) * cy + kx2 * sx * sy)) / (kx2 ** 2 - 4 * kx2 * ky2)
        T[2, 2, 5] = (ky2 * L * sy) / 2
        T[2, 3, 5] = I34 * ky2
        T[3, 0, 2] = (2 * cy * K2 * (kx2 - 2 * ky2) * sx - 2 * (1 + cx) * K2 * ky2 * sy) / (kx2 - 4 * ky2)
        T[3, 0, 3] = (2 * K2 * ((-1 + cx) * cy + (kx2 - 2 * ky2) * sx * sy)) / (kx2 - 4 * ky2)
        T[3, 1, 2] = (-2 * K2 * ((-1 + cx) * cy * (kx2 - 2 * ky2) + kx2 * ky2 * sx * sy)) / (kx2 ** 2 - 4 * kx2 * ky2)
        T[3, 1, 3] = (2 * K2 * (cy * kx2 * sx - (cx * (kx2 - 2 * ky2) + 2 * ky2) * sy)) / (kx2 ** 2 - 4 * kx2 * ky2)
        T[3, 2, 5] = (ky2 * (cy * L + sy)) / 2
        T[3, 3, 5] = (ky2 * L * sy) / 2

    if kx2 == 0 and ky2 == 0:
        T[0, 0, 0] = -(K2 * L ** 2) / 2
        T[0, 0, 1] = -(K2 * L ** 3) / 3
        T[0, 0, 5] = 0
        T[0, 1, 1] = -(K2 * L ** 4) / 12
        T[0, 1, 5] = 0
        T[0, 5, 5] = 0
        T[0, 2, 2] = (K2 * L ** 2) / 2
        T[0, 2, 3] = (K2 * L ** 3) / 3
        T[0, 3, 3] = (K2 * L ** 4) / 12
        T[1, 0, 0] = -(K2 * L)
        T[1, 0, 1] = -(K2 * L ** 2)
        T[1, 0, 5] = 0
        T[1, 1, 1] = -(K2 * L ** 3) / 3
        T[1, 1, 5] = 0
        T[1, 5, 5] = 0
        T[1, 2, 2] = K2 * L
        T[1, 2, 3] = K2 * L ** 2
        T[1, 3, 3] = (K2 * L ** 3) / 3
        T[2, 0, 2] = K2 * L ** 2
        T[2, 0, 3] = (-2 * K2 * L ** 3) / 3
        T[2, 1, 2] = (K2 * L ** 3) / 3
        T[2, 1, 3] = (K2 * L ** 4) / 6
        T[2, 2, 5] = 0
        T[2, 3, 5] = 0
        T[3, 0, 2] = 2 * K2 * L
        T[3, 0, 3] = 2 * K2 * L ** 2
        T[3, 1, 2] = 0
        T[3, 1, 3] = -(K2 * L ** 3) / 3
        T[3, 2, 5] = 0
        T[3, 3, 5] = 0

    if kx2 == 4 * ky2 and h != 0:
        T[0, 0, 0] = ((-1 + cy ** 2) * (h ** 3 + K2 + 2 * h * (K1 - 2 * ky2)) - (h ** 3 + 2 * h * K1 + K2 + 4 * (
                    h ** 3 + h * (-1 + 2 * K1) + K2) * ky2 + 8 * h * ky2 ** 2) * sy ** 2) / (12 * ky2)
        T[0, 0, 1] = (h * sy * (
                    6 + 2 * (-1 + cy ** 2 - sy ** 2) + ((h ** 3 + 2 * h * K1 + K2) * (-1 + cy ** 2 - sy ** 2)) / (
                        h * ky2))) / 6
        T[0, 0, 5] = (-((-1 + cy ** 2) * h * (h ** 3 + K2 + 2 * h * (K1 - 2 * ky2))) - 6 * ky2 * (
                    h ** 4 + h * K2 + 2 * h ** 2 * (K1 - 2 * ky2) - 2 * K1 * ky2) * L * sy + h * (
                                  h ** 3 + 2 * h * K1 + K2 + 4 * (
                                      h ** 3 + h * (-1 + 2 * K1) + K2) * ky2 + 8 * h * ky2 ** 2) * sy ** 2) / (
                                 24 * ky2 ** 2)
        T[0, 1, 1] = ((-1 + cy ** 2) * (h ** 3 + 2 * h * K1 + K2 - h * ky2) + (
                    h ** 3 * (-1 + 2 * ky2) + K2 * (-1 + 2 * ky2) + h * (
                        -2 * K1 + ky2 + 4 * K1 * ky2 + 4 * ky2 ** 2)) * sy ** 2) / (24 * ky2 ** 2)
        T[0, 1, 5] = (4 * h ** 2 * ky2 * sy * (1 - cy ** 2 + sy ** 2) + 6 * (2 * h ** 2 + K1) * ky2 * (
                    -(cy ** 2 * L) + sy + L * sy ** 2) + h * (h ** 3 + 2 * h * K1 + K2) * (
                                  -sy + (3 * L - 2 * sy) * (cy - sy) * (cy + sy))) / (48 * ky2 ** 2)
        T[0, 5, 5] = (h * ((-1 + cy ** 2) * (
                    h ** 4 + h * K2 + h ** 2 * (2 * K1 - 7 * ky2) + 3 * ky2 * (-K1 + 4 * ky2)) + 3 * ky2 * (
                                       h ** 4 + h * K2 + 2 * h ** 2 * (K1 - 2 * ky2) - 2 * K1 * ky2) * L * sy - (
                                       h * (h ** 3 + 2 * h * K1 + K2) - 3 * K1 * ky2 + h * (
                                           h ** 3 + h * (-7 + 2 * K1) + K2) * ky2 + 2 * (
                                                   6 + h ** 2) * ky2 ** 2) * sy ** 2)) / (48 * ky2 ** 3)
        T[0, 2, 2] = ((h ** 5 - 5 * h ** 2 * K2 - 5 * h ** 3 * ky2 + 50 * K2 * ky2) * (
                    -1 + cy ** 2 - sy ** 2) + 50 * K2 * ky2 * sin((h * L) / sqrt(5)) ** 2) / (
                                 40 * (h ** 2 - 5 * ky2) * ky2)
        T[0, 2, 3] = (5 * K2 * (2 * h * sy - sqrt(5) * sin((2 * h * L) / sqrt(5)))) / (4 * (h ** 3 - 5 * h * ky2))
        T[0, 3, 3] = (h ** 2 * (h ** 3 - 5 * K2 - 5 * h * ky2) * (-1 + cy ** 2 - sy ** 2) - 50 * K2 * ky2 * sin(
            (h * L) / sqrt(5)) ** 2) / (8 * h ** 2 * (h ** 2 - 5 * ky2) * ky2)
        T[1, 0, 0] = (h * sy * (1 - cy ** 2 + sy ** 2 + 12 * ky2 * (cy - sy) * (cy + sy) - (
                    (h ** 3 + 2 * h * K1 + K2) * (1 + 2 * cy ** 2 - 2 * sy ** 2)) / h)) / 3
        T[1, 0, 1] = (-((-1 + cy ** 2) * (h ** 3 + K2 + 2 * h * (K1 + ky2 + 3 * cy ** 2 * ky2))) + (
                    h ** 3 + 2 * h * K1 + K2 - 4 * (
                        h - 3 * cy ** 2 * h + 2 * h ** 3 + 4 * h * K1 + 2 * K2) * ky2 + 8 * h * ky2 ** 2) * sy ** 2 - 6 * h * ky2 * sy ** 4) / (
                                 6 * ky2)
        T[1, 0, 5] = (-3 * cy ** 2 * (h ** 4 + h * K2 + 2 * h ** 2 * (K1 - 2 * ky2) - 2 * K1 * ky2) * L + (
                    (-1 + 4 * cy ** 2) * h * (h ** 3 + 2 * h * K1 + K2) - 16 * (
                        -1 + cy ** 2) * h ** 2 * ky2 + 6 * K1 * ky2) * sy + 3 * (h ** 4 + h * K2 + 2 * h ** 2 * (
                    K1 - 2 * ky2) - 2 * K1 * ky2) * L * sy ** 2 - 4 * h * (
                                  h ** 3 + K2 + 2 * h * (K1 - 2 * ky2)) * sy ** 3) / (12 * ky2)
        T[1, 1, 1] = (sy * (h ** 3 * (-1 + cy ** 2 - sy ** 2) + K2 * (-1 + cy ** 2 - sy ** 2) + h * (
                    ky2 - 4 * cy ** 2 * ky2 + 4 * ky2 * sy ** 2 + 2 * K1 * (-1 + cy ** 2 - sy ** 2)))) / (6 * ky2)
        T[1, 1, 5] = ((-1 + cy ** 2) * h * (h ** 3 + K2 + 2 * h * (K1 + ky2 + 3 * cy ** 2 * ky2)) - 6 * ky2 * (
                    h ** 4 + h * K2 + 2 * h ** 2 * (K1 - 2 * ky2) - 2 * K1 * ky2) * L * sy + h * (
                                  h ** 3 * (-1 + 8 * ky2) + K2 * (-1 + 8 * ky2) - 2 * h * (
                                      K1 - 8 * K1 * ky2 + 2 * ky2 * (
                                          -1 + 3 * cy ** 2 + 2 * ky2))) * sy ** 2 + 6 * h ** 2 * ky2 * sy ** 4) / (
                                 24 * ky2 ** 2)
        T[1, 5, 5] = (h * (-48 * ky2 ** 2 * sy - 8 * h ** 2 * ky2 * sy * (1 - cy ** 2 + sy ** 2) + 6 * (
                    2 * h ** 2 + K1) * ky2 * (-(cy ** 2 * L) + sy + L * sy ** 2) + h * (h ** 3 + 2 * h * K1 + K2) * (
                                       -sy + (3 * L - 2 * sy) * (cy - sy) * (cy + sy)))) / (48 * ky2 ** 2)
        T[1, 2, 2] = (-2 * (h ** 5 - 5 * h ** 2 * K2 - 5 * h ** 3 * ky2 + 50 * K2 * ky2) * sy + 5 * sqrt(
            5) * h * K2 * sin((2 * h * L) / sqrt(5))) / (20 * (h ** 2 - 5 * ky2))
        T[1, 2, 3] = (-5 * K2 * (-cy ** 2 + sy ** 2 + cos((2 * h * L) / sqrt(5)))) / (2 * (h ** 2 - 5 * ky2))
        T[1, 3, 3] = -(2 * h * (h ** 3 - 5 * K2 - 5 * h * ky2) * sy + 5 * sqrt(5) * K2 * sin((2 * h * L) / sqrt(5))) / (
                    4 * (h ** 3 - 5 * h * ky2))
        T[2, 0, 2] = ((h ** 3 + 5 * h * K1 + 5 * K2) * (-1 + cy ** 2 - sy ** 2) * cos((h * L) / sqrt(5)) + 2 * sqrt(
            5) * h * (K2 + h * (K1 + ky2)) * sy * sin((h * L) / sqrt(5))) / (2 * (h ** 2 - 5 * ky2))
        T[2, 0, 3] = (-10 * h * (K2 + h * (K1 + ky2)) * sy * cos((h * L) / sqrt(5)) + sqrt(5) * (
                    h ** 3 + 5 * h * K1 + 5 * K2) * (1 + cy ** 2 - sy ** 2) * sin((h * L) / sqrt(5))) / (
                                 2 * (h ** 3 - 5 * h * ky2))
        T[2, 1, 2] = (2 * h * (h ** 3 + 5 * h * K1 + 5 * K2) * ky2 * sy * cos((h * L) / sqrt(5)) + sqrt(5) * (
                    -10 * h * K1 * ky2 - 10 * K2 * ky2 + h ** 2 * K2 * (1 - cy ** 2 + sy ** 2) + h ** 3 * (
                        K1 - cy ** 2 * K1 - ky2 - cy ** 2 * ky2 + (K1 + ky2) * sy ** 2)) * sin((h * L) / sqrt(5))) / (
                                 4 * h * (h ** 2 - 5 * ky2) * ky2)
        T[2, 1, 3] = (5 * h * (K2 + h * (K1 + ky2)) * (-1 + cy ** 2 - sy ** 2) * cos((h * L) / sqrt(5)) + 2 * sqrt(
            5) * (h ** 3 + 5 * h * K1 + 5 * K2) * ky2 * sy * sin((h * L) / sqrt(5))) / (
                                 4 * h * (h ** 2 - 5 * ky2) * ky2)
        T[2, 2, 5] = (-5 * h * (h ** 3 + 5 * h * K1 + 5 * K2) * (-1 + cy ** 2 - sy ** 2) * cos(
            (h * L) / sqrt(5)) + 2 * sqrt(5) * (
                                  (h ** 2 - 5 * ky2) * (5 * h * K1 + 5 * K2 + 2 * h * ky2) * L - 5 * h ** 2 * (
                                      K2 + h * (K1 + ky2)) * sy) * sin((h * L) / sqrt(5))) / (
                                 40 * (h ** 2 - 5 * ky2) * ky2)
        T[2, 3, 5] = (4 * h * I34 * (h ** 2 - 5 * ky2) * (5 * h * K1 + 5 * K2 + 2 * h * ky2) + 50 * h * (
                    K2 + h * (K1 + ky2)) * sy * cos((h * L) / sqrt(5)) - 5 * sqrt(5) * (
                                  h ** 3 * (-1 + cy ** 2 - sy ** 2) + 5 * K2 * (1 + cy ** 2 - sy ** 2) + 5 * h * (
                                      K1 + cy ** 2 * K1 + 2 * ky2 - K1 * sy ** 2)) * sin((h * L) / sqrt(5))) / (
                                 40 * (h ** 2 - 5 * ky2) * ky2)
        T[3, 0, 2] = (10 * (h ** 2 * (h * K1 + K2) - (h ** 3 + 10 * h * K1 + 10 * K2) * ky2) * sy * cos(
            (h * L) / sqrt(5)) + sqrt(5) * h * (h ** 3 + 5 * h * K1 + 5 * K2) * (1 + cy ** 2 - sy ** 2) * sin(
            (h * L) / sqrt(5))) / (10 * (h ** 2 - 5 * ky2))
        T[3, 0, 3] = (-(
                    h * (h ** 3 + 5 * h * K1 + 5 * K2) * (-1 + cy ** 2 - sy ** 2) * cos((h * L) / sqrt(5))) + 2 * sqrt(
            5) * (h ** 2 * (h * K1 + K2) - (h ** 3 + 10 * h * K1 + 10 * K2) * ky2) * sy * sin((h * L) / sqrt(5))) / (
                                 2 * (h ** 3 - 5 * h * ky2))
        T[3, 1, 2] = (-5 * (h ** 2 * (h * K1 + K2) - (h ** 3 + 10 * h * K1 + 10 * K2) * ky2) * (
                    -1 + cy ** 2 - sy ** 2) * cos((h * L) / sqrt(5)) + 2 * sqrt(5) * h * (
                                  h ** 3 + 5 * h * K1 + 5 * K2) * ky2 * sy * sin((h * L) / sqrt(5))) / (
                                 20 * (h ** 2 - 5 * ky2) * ky2)
        T[3, 1, 3] = (-2 * h * (h ** 3 + 5 * h * K1 + 5 * K2) * ky2 * sy * cos((h * L) / sqrt(5)) + sqrt(5) * (
                    10 * h * K1 * ky2 * (cy - sy) * (cy + sy) + 10 * K2 * ky2 * (cy - sy) * (cy + sy) + h ** 2 * K2 * (
                        1 - cy ** 2 + sy ** 2) + h ** 3 * (
                                K1 - cy ** 2 * K1 + ky2 + cy ** 2 * ky2 + (K1 - ky2) * sy ** 2)) * sin(
            (h * L) / sqrt(5))) / (4 * h * (h ** 2 - 5 * ky2) * ky2)
        T[3, 2, 5] = (2 * h * ((h ** 2 - 5 * ky2) * (5 * h * K1 + 5 * K2 + 2 * h * ky2) * L - 5 * h ** 2 * (
                    h * K1 + K2) * sy + 5 * (h ** 3 + 10 * h * K1 + 10 * K2) * ky2 * sy) * cos(
            (h * L) / sqrt(5)) + sqrt(5) * (-50 * K2 * ky2 - 10 * h * ky2 * (5 * K1 + 2 * ky2) + h ** 5 * (
                    1 - cy ** 2 + sy ** 2) + 5 * h ** 2 * K2 * (1 - cy ** 2 + sy ** 2) + h ** 3 * (
                                                        -6 * ky2 + 5 * K1 * (1 - cy ** 2 + sy ** 2))) * sin(
            (h * L) / sqrt(5))) / (40 * (h ** 2 - 5 * ky2) * ky2)
        T[3, 3, 5] = (5 * h * (h ** 3 + 5 * h * K1 + 5 * K2) * (-1 + cy ** 2 - sy ** 2) * cos(
            (h * L) / sqrt(5)) + 2 * sqrt(5) * (
                                  (h ** 2 - 5 * ky2) * (5 * h * K1 + 5 * K2 + 2 * h * ky2) * L - 5 * h ** 2 * (
                                      h * K1 + K2) * sy + 5 * (h ** 3 + 10 * h * K1 + 10 * K2) * ky2 * sy) * sin(
            (h * L) / sqrt(5))) / (40 * (h ** 2 - 5 * ky2) * ky2)

    return T

