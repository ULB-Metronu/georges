from __future__ import annotations
from typing import Optional
import numpy as _np
from .input import Input as _Input
from .beam import Beam as _Beam
from .observers import Observer


def track(beamline: _Input,
          beam: _Beam,
          observer: Optional[Observer] = None,
          check_apertures: bool = True
          ):
    """

    Args:
        beamline:
        beam:
        observer:
        check_apertures:

    Returns:

    """
    global_parameters = [
        beam.kinematics.beta
    ]
    b1 = beam.distribution
    b2 = _np.zeros(b1.shape)
    for e in beamline.sequence:
        b1, b2 = e.propagate(b1, b2, global_parameters)
        if check_apertures:
            b2, b1 = e.check_aperture(b1, b2)
        else:
            b2, b1 = b1, b2
        if observer is not None:
            observer(e, b1, b2)


def twiss(beamline: _Input, periodic: bool = False):
    """

    Args:
        beamline:
        periodic:

    Returns:

    """
    pass
