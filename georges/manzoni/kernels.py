"""
TODO
"""
import numpy as _np
import numba as _nb
from numba import njit


@njit(parallel=True, nogil=True)
def batched_vector_matrix(b1: _np.ndarray, b2: _np.ndarray, matrix: _np.ndarray):
    """

    Args:
        b1: a numpy array containing all the particles
        b2: explicit destination for the result
        matrix: the transfer matrix as a numpy array

    Returns:
        the destination (result) array
    """
    return b1, _np.dot(b1, matrix.T, out=b2)


@njit(parallel=True, nogil=True)
def batched_vector_tensor(b1: _np.ndarray, b2: _np.ndarray, tensor: _np.ndarray):
    """

    Args:
        b1:
        b2:
        tensor:

    Returns:

    """
    for l in _nb.prange(b1.shape[0]):
        for i in range(tensor.shape[0]):
            s = 0.0
            for j in range(tensor.shape[1]):
                for k in range(tensor.shape[2]):
                    s += tensor[i, j, k] * b1[l, j] * b1[l, k]
            b2[l, i] = s
    return b1, b2


@njit(parallel=True, nogil=True)
def batched_vector_matrix_tensor(b1: _np.ndarray, b2: _np.ndarray, matrix: _np.ndarray, tensor: _np.ndarray):
    """

    Args:
        b1:
        b2:
        matrix:
        tensor:

    Returns:

    """
    for l in range(b1.shape[0]):
        for i in range(tensor.shape[0]):
            s = 0
            for j in range(tensor.shape[1]):
                s += matrix[i, j] * b1[l, j]
                for k in range(tensor.shape[2]):
                    s += tensor[i, j, k] * b1[l, j] * b1[l, k]
            b2[l, i] = s
    return b1, b2


@njit
def matrix_matrix(m1, m2):
    return _np.matmul(m1, m2)
