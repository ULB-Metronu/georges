"""
TODO
"""
import numpy as _np
from numba import njit


@njit
def circular_aperture_check(b1, kargs: _np.ndarray):
    """

    Args:
        b1:
        kargs:

    Returns:

    """
    return (b1[:, 0] ** 2 + b1[:, 2] ** 2) < kargs[0] ** 2


@njit
def elliptical_aperture_check(b1, kargs: _np.ndarray):
    """

    Args:
        b1:
        kargs:

    Returns:

    """
    return (b1[:, 0] / kargs[0]) ** 2 + (b1[:, 2] / kargs[1] ** 2) < 1


@njit
def rectangular_aperture_check(b1, kargs: _np.ndarray):
    """

    Args:
        b1:
        kargs:

    Returns:

    """
    return _np.multiply(_np.abs(b1[:, 0]) < kargs[0], _np.abs(b1[:, 2]) < kargs[1])


@njit
def phase_space_aperture_check(b1, kargs: _np.ndarray):
    """

    Args:
        b1:
        kargs:

    Returns:

    """
    return ((b1[:, 0] / kargs[0]) ** 2 + (b1[:, 1] / kargs[1]) ** 2) < 1 and (
        (b1[:, 2] / kargs[2]) ** 2 + (b1[:, 3] / kargs[3]) ** 2
    ) < 1
