"""
TODO
"""
import numpy as _np
from numba import njit


@njit(parallel=True)
def circular_aperture_check(b1, kargs: _np.ndarray):
    """

    Args:
        b1:
        kargs:

    Returns:

    """
    return (b1[:, 0] ** 2 + b1[:, 2] ** 2) < kargs[0] ** 2


@njit(parallel=True)
def elliptical_aperture_check(b1, kargs: _np.ndarray):
    """

    Args:
        b1:
        kargs:

    Returns:

    """
    return (b1[:, 0] ** 2 + b1[:, 2] ** 2) < kargs[0] ** 2


@njit(parallel=True)
def rectangular_aperture_check(b1, kargs: _np.ndarray):
    """

    Args:
        b1:
        kargs:

    Returns:

    """
    return (b1[:, 0] ** 2 + b1[:, 2] ** 2) < kargs[0] ** 2


@njit(parallel=True)
def phase_space_aperture_check(b1, kargs: _np.ndarray):
    """

    Args:
        b1:
        kargs:

    Returns:

    """
    return (b1[:, 0] ** 2 + b1[:, 2] ** 2) < kargs[0] ** 2
