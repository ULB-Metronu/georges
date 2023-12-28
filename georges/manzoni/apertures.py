"""
The file `apertures.py` is for the `collimator` elements and contains the implementation of the
selection of the particle according to the type of collimator (circular, elliptical, rectangular or phase-space).
"""
import numpy as _np
from numba import njit


@njit
def circular_aperture_check(b1, kargs: _np.ndarray):
    """
    The aperture is circular and defined by its radius.

    Args:
        b1: The beam before the collimator
        kargs: Radius of the aperture.

    Returns:
        The collimated beam

    """
    return (b1[:, 0] ** 2 + b1[:, 2] ** 2) < kargs[0] ** 2


@njit
def elliptical_aperture_check(b1, kargs: _np.ndarray):
    """
    The aperture is elliptical and defined by the semi-axes of the elliptical aperture.

    Args:
        b1: The beam before the collimator
        kargs: Semi-axes of the aperture.

    Returns:
        The collimated beam

    """
    return (b1[:, 0] / kargs[0]) ** 2 + (b1[:, 2] / kargs[1] ** 2) < 1


@njit
def rectangular_aperture_check(b1, kargs: _np.ndarray):
    """
    The aperture is rectangular and defined by its horizontal and vertical half apertures.

    Args:
        b1: The beam before the collimator
        kargs: horizontal and vertical half apertures

    Returns:
        The collimated beam

    """
    return _np.multiply(_np.abs(b1[:, 0]) < kargs[0], _np.abs(b1[:, 2]) < kargs[1])


@njit
def phase_space_aperture_check(b1, kargs: _np.ndarray):
    """
    The phase space aperture allows the user to cut either the position or the momentum of the beam.

    Args:
        b1: The beam before the collimator
        kargs: Radius of the aperture in each coordinates system (position and momentum).

    Returns:
        The collimated beam

    """
    return ((b1[:, 0] / kargs[0]) ** 2 + (b1[:, 1] / kargs[1]) ** 2) < 1 and (
        (b1[:, 2] / kargs[2]) ** 2 + (b1[:, 3] / kargs[3]) ** 2
    ) < 1
