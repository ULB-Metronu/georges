from numba import njit
import numpy as np
from numpy import cos, sin, cosh, sinh, sqrt


@njit
def compute_mad_quadrupole_matrix(L: float, K1: float, beta: float) -> np.ndarray:
    gamma = 1 / sqrt(1-beta**2)
    R = np.zeros((6, 6))

    # Setup basic notations - Horizontal plane
    if K1 == 0:
        cx = 1
        sx = L
    elif K1 > 0:
        cx = cos(sqrt(K1)*L)
        sx = sin(sqrt(K1)*L)/sqrt(K1)
    else:
        cx = cosh(sqrt(-K1)*L)
        sx = sinh(sqrt(-K1)*L)/sqrt(-K1)
    
    # Setup basic notations - Vertical plane
    if K1 == 0:
        cy = 1
        sy = L
    elif K1 > 0:
        cy = cosh(sqrt(K1)*L)
        sy = sinh(sqrt(K1)*L)/sqrt(K1)
    else:
        cy = cos(sqrt(-K1)*L)
        sy = sin(sqrt(-K1)*L)/sqrt(-K1)
    
    # Definition of the matrix elements
    R[0, 0] = cx
    R[0, 1] = sx
    R[1, 0] = -(K1*sx)
    R[1, 1] = cx
    R[2, 2] = cy
    R[2, 3] = sy
    R[3, 2] = K1*sy
    R[3, 3] = cy
    R[4, 4] = 1
    R[5, 5] = 1
    
    return R


@njit
def compute_mad_quadrupole_tensor(L: float, K1: float, beta: float) -> np.ndarray:
    gamma = 1 / sqrt(1-beta**2)
    T = np.zeros((6, 6, 6))

    # Setup basic notations - Horizontal plane
    if K1 == 0:
        cx = 1
        sx = L
    elif K1 > 0:
        cx = cos(sqrt(K1)*L)
        sx = sin(sqrt(K1)*L)/sqrt(K1)
    else:
        cx = cosh(sqrt(-K1)*L)
        sx = sinh(sqrt(-K1)*L)/sqrt(-K1)
    
    # Setup basic notations - Vertical plane
    if K1 == 0:
        cy = 1
        sy = L
    elif K1 > 0:
        cy = cosh(sqrt(K1)*L)
        sy = sinh(sqrt(K1)*L)/sqrt(K1)
    else:
        cy = cos(sqrt(-K1)*L)
        sy = sin(sqrt(-K1)*L)/sqrt(-K1)
    
    # Definition of the tensor elements
    T[0, 0, 5] = (K1*L*sx)/(2*beta)
    T[0, 1, 5] = (-(cx*L) - sx)/(2*beta)
    T[1, 0, 5] = -(K1*(-(cx*L) + sx))/(2*beta)
    T[1, 1, 5] = (K1*L*sx)/(2*beta)
    T[2, 2, 5] = -(K1*L*sy)/(2*beta)
    T[2, 3, 5] = (-(cy*L) - sy)/(2*beta)
    T[3, 2, 3] = (K1*(-(cy*L) + sy))/(2*beta)
    T[3, 3, 5] = -(K1*L*sy)/(2*beta)
    T[4, 0, 0] = -(K1*(L - cx*sx))/(4*beta)
    T[4, 0, 1] = (K1*sx**2)/(2*beta)
    T[4, 1, 1] = -(L + cx*sx)/(4*beta)
    T[4, 2, 2] = (K1*(L - cy*sy))/(4*beta)
    T[4, 2, 3] = (K1*sy**2)/(2*beta)
    T[4, 3, 3] = -(L + cy*sy)/(4*beta)
    T[4, 5, 5] = (-3*L)/(2*beta**3*gamma**2)
    
    return T

