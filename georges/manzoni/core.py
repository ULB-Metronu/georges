from __future__ import annotations
from typing import Optional
import numpy as _np
from .input import Input as _Input
from .observers import Observer


def track(beamline: _Input, beam: _np.ndarray, observer: Optional[Observer] = None):
    """

    Args:
        beamline:
        beam:
        observer:

    Returns:

    """
    b1 = beam
    b2 = _np.zeros(beam.shape)
    for e in beamline.sequence:
        b1, b2 = e.propagate(b1, b2)
        #b2, b1 = e.check_aperture(b1, b2)
        #b2 = _np.zeros(b1.shape)
        if observer is not None:
            observer(e, b1, b2)
        b1 = b2
        b2 = _np.zeros(b1.shape)


def twiss(beamline: _Input, periodic: bool = False):
    """

    Args:
        beamline:
        periodic:

    Returns:

    """
    pass
