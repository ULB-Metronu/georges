"""
TODO
"""
import numpy as _np
import numba as _nb
from numba import njit


@njit(parallel=True)
def batched_vector_matrix(b1: _np.ndarray, b2: _np.ndarray, kargs: _np.ndarray):
    """

    Args:
        b1: a numpy array containing all the particles
        b2: explicit destination for the result
        kargs: the transfer matrix as a numpy array

    Returns:
        the destination (result) array
    """
    return _np.dot(b1, kargs.T, out=b2)


@njit
def matrix_matrix(m1, m2):
    return _np.matmul(m1, m2)


def drift_kick_drift():
    pass


def matrix_kick_matrix():
    pass


@njit(parallel=True)
def drift2(b1: _np.ndarray, b2: _np.ndarray, length: float):
    for i in _nb.prange(0, b1.shape[0]):
        b2[i, 0] = length * b1[i, 1]


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
