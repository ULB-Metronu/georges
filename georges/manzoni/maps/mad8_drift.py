import numba as _nb
import numpy as _np
from numba import njit


@njit(parallel=True)
def drift2(b1: _np.ndarray, b2: _np.ndarray, length: float):
    for i in _nb.prange(0, b1.shape[0]):
        b2[i, 0] = length * b1[i, 1]


@njit(parallel=True)
def drift4(b1: _np.ndarray, b2: _np.ndarray, kargs: _np.ndarray):
    """
    Performance tests on MBP i9 show that the explicit loop is faster, both at compile-time and at runtime
    than the nb broadcast variant.

    Args:
        b1: a numpy array containing all the particles
        b2: explicit destination for the result
        kargs: the drift length (in meters)

    Returns:
        the destination (result) array
    """
    for i in _nb.prange(0, b1.shape[0]):
        b2[i, 0] = b1[i, 0] + kargs[0] * b1[i, 1]
        b2[i, 1] = b1[i, 1]
        b2[i, 2] = b1[i, 2] + kargs[0] * b1[i, 3]
        b2[i, 3] = b1[i, 3]
    return b2


@njit(parallel=True)
def drift5(b1: _np.ndarray, b2: _np.ndarray, kargs: _np.ndarray):
    """
    Performance tests on MBP i9 show that the explicit loop is faster, both at compile-time and at runtime
    than the nb broadcast variant.

    Args:
        b1: a numpy array containing all the particles
        b2: explicit destination for the result
        kargs: the drift length (in meters)

    Returns:
        the destination (result) array
    """
    for i in _nb.prange(0, b1.shape[0]):
        b2[i, 0] = b1[i, 0] + kargs[0] * b1[i, 1]
        b2[i, 1] = b1[i, 1]
        b2[i, 2] = b1[i, 2] + kargs[0] * b1[i, 3]
        b2[i, 3] = b1[i, 3]
        b2[i, 4] = b1[i, 4]
    return b2


@njit(parallel=True)
def drift6(b1: _np.ndarray, b2: _np.ndarray, length: float):
    """

    .. note:: Performance tests on MBP i9 show that the explicit loop is faster, both at compile-time and at runtime
    than the nb broadcast variant.

    Args:
        b1: a numpy array containing all the particles
        b2: explicit destination for the result
        length: the drift length (in meters)

    Returns:
        the destination (result) array
    """
    for i in _nb.prange(0, b1.shape[0]):
        b2[i, 0] = b1[i, 0] + length * b1[i, 1]
        b2[i, 1] = b1[i, 1]
        b2[i, 2] = b1[i, 2] + length * b1[i, 3]
        b2[i, 3] = b1[i, 3]
        b2[i, 4] = b1[i, 4]
        b2[i, 5] = b1[i, 5]
    return b1, b2


@njit(parallel=True, fastmath=True)
def compute_mad_drift_matrix(element_parameters: list, **_) -> _np.ndarray:
    L: float = element_parameters[0]
    R = _np.zeros((6, 6))

    # Definition of the matrix elements
    R[0, 0] = 1
    R[0, 1] = L
    R[1, 1] = 1
    R[2, 2] = 1
    R[2, 3] = L
    R[3, 3] = 1
    R[4, 4] = 1
    R[5, 5] = 1

    return R


@njit(parallel=True, fastmath=True)
def compute_mad_drift_tensor(element_parameters: list, **_) -> _np.ndarray:
    T = _np.zeros((6, 6, 6))

    return T
