from numba import njit
import numpy as np


@njit(cache=True)
def compute_mad_quadrupole_matrix(element_parameters: list, **_) -> np.ndarray:
    L: float = element_parameters[0]
    k1: float = element_parameters[1]
    R = np.zeros((6, 6))

    # Setup basic notations - Horizontal plane
    if k1 == 0:
        cx = 1
        sx = L
    elif k1 > 0:
        cx = np.cos(np.sqrt(k1)*L)
        sx = np.sin(np.sqrt(k1)*L)/np.sqrt(k1)
    else:
        cx = np.cosh(np.sqrt(-k1)*L)
        sx = np.sinh(np.sqrt(-k1)*L)/np.sqrt(-k1)
    
    # Setup basic notations - Vertical plane
    if k1 == 0:
        cy = 1
        sy = L
    elif k1 > 0:
        cy = np.cosh(np.sqrt(k1)*L)
        sy = np.sinh(np.sqrt(k1)*L)/np.sqrt(k1)
    else:
        cy = np.cos(np.sqrt(-k1)*L)
        sy = np.sin(np.sqrt(-k1)*L)/np.sqrt(-k1)
    
    # Definition of the matrix elements
    R[0, 0] = cx
    R[0, 1] = sx
    R[1, 0] = -(k1*sx)
    R[1, 1] = cx
    R[2, 2] = cy
    R[2, 3] = sy
    R[3, 2] = k1*sy
    R[3, 3] = cy
    R[4, 4] = 1
    R[5, 5] = 1
    
    return R


@njit(cache=True)
def compute_mad_quadrupole_tensor(element_parameters: list, global_parameters: list) -> np.ndarray:
    L: float = element_parameters[0]
    k1: float = element_parameters[1]
    beta: float = global_parameters[0]
    gamma = 1 / np.sqrt(1-beta**2)
    T = np.zeros((6, 6, 6))

    # Setup basic notations - Horizontal plane
    if k1 == 0:
        cx = 1
        sx = L
    elif k1 > 0:
        cx = np.cos(np.sqrt(k1)*L)
        sx = np.sin(np.sqrt(k1)*L)/np.sqrt(k1)
    else:
        cx = np.cosh(np.sqrt(-k1)*L)
        sx = np.sinh(np.sqrt(-k1)*L)/np.sqrt(-k1)
    
    # Setup basic notations - Vertical plane
    if k1 == 0:
        cy = 1
        sy = L
    elif k1 > 0:
        cy = np.cosh(np.sqrt(k1)*L)
        sy = np.sinh(np.sqrt(k1)*L)/np.sqrt(k1)
    else:
        cy = np.cos(np.sqrt(-k1)*L)
        sy = np.sin(np.sqrt(-k1)*L)/np.sqrt(-k1)
    
    # Definition of the tensor elements
    T[0, 0, 5] = (k1*L*sx)/(2*beta)
    T[0, 1, 5] = (-(cx*L) - sx)/(2*beta)
    T[1, 0, 5] = -(k1*(-(cx*L) + sx))/(2*beta)
    T[1, 1, 5] = (k1*L*sx)/(2*beta)
    T[2, 2, 5] = -(k1*L*sy)/(2*beta)
    T[2, 3, 5] = (-(cy*L) - sy)/(2*beta)
    T[3, 2, 5] = (k1*(-(cy*L) + sy))/(2*beta)
    T[3, 3, 5] = -(k1*L*sy)/(2*beta)
    T[4, 0, 0] = -(k1*(L - cx*sx))/(4*beta)
    T[4, 0, 1] = (k1*sx**2)/(2*beta)
    T[4, 1, 1] = -(L + cx*sx)/(4*beta)
    T[4, 2, 2] = (k1*(L - cy*sy))/(4*beta)
    T[4, 2, 3] = (k1*sy**2)/(2*beta)
    T[4, 3, 3] = -(L + cy*sy)/(4*beta)
    T[4, 5, 5] = (-3*L)/(2*beta**3*gamma**2)
    
    return T
