"""
The file `kernels.py` contains the loops that are the core of the particles propagation based on
their coordinates. Different batches are available, to allow a matrix (order 1) propagation,
a tensor (order 2) propagation or a matrix followed by a tensor (orders 1+2) propagations.
"""
import numba as _nb
import numpy as _np
from numba import njit


@njit(nogil=True)  # no parallel jit found
def batched_vector_matrix(b1: _np.ndarray, b2: _np.ndarray, matrix: _np.ndarray):  # pragma: no cover
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
def batched_vector_tensor(b1: _np.ndarray, b2: _np.ndarray, tensor: _np.ndarray):  # pragma: no cover
    """

    Args:
        b1:
        b2:
        tensor:

    Returns:

    """
    for h in _nb.prange(b1.shape[0]):
        for i in range(tensor.shape[0]):
            s = 0.0
            for j in range(tensor.shape[1]):
                for k in range(j, tensor.shape[2]):  # Assume upper triangular matrix Mjk = Ti::
                    s += tensor[i, j, k] * b1[h, j] * b1[h, k]
            b2[h, i] = s
    return b1, b2


@njit(nogil=True)
def batched_vector_matrix_tensor(
    b1: _np.ndarray,
    b2: _np.ndarray,
    matrix: _np.ndarray,
    tensor: _np.ndarray,
):  # pragma: no cover
    """

    Args:
        b1:
        b2:
        matrix:
        tensor:

    Returns:

    """
    for h in range(b1.shape[0]):
        for i in range(tensor.shape[0]):
            s = 0
            for j in range(tensor.shape[1]):
                s += matrix[i, j] * b1[h, j]
                for k in range(j, tensor.shape[2]):  # Assume upper triangular matrix Mjk = Ti::
                    s += tensor[i, j, k] * b1[h, j] * b1[h, k]
            b2[h, i] = s
    return b1, b2


@njit
def matrix_matrix(m1, m2):  # pragma: no cover
    return _np.matmul(m1, m2)
