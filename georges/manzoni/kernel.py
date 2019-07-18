import numpy as np
import numba as nb
from numba import njit


@njit(parallel=True)
def batched_vector_matrix(p, m):
    return np.dot(p, m.T)


@njit
def matrix_matrix(m1, m2):
    return np.matmul(m1, m2)

def drift_kick_drift():
    pass

def matrix_kick_matrix():
    pass
