from __future__ import annotations
from typing import TYPE_CHECKING, Optional
import numpy as _np
if TYPE_CHECKING:
    from .input import Input as _Input
    from .beam import Beam as _Beam
    from .observers import Observer as _Observer


def track(beamline: _Input,
          beam: _Beam,
          observer: Optional[_Observer] = None,
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
    b1 = _np.copy(beam.distribution)
    b2 = _np.zeros(b1.shape)
    for e in beamline.sequence:
        if check_apertures_entry:
            b1, b2 = e.check_aperture(b1, b2)
            if b1.shape != b2.shape:
                b1 = _np.zeros(b2.shape)
            if b2.shape[0] == 0:
                break
        if check_apertures_exit:
            if b1.shape != b2.shape:
                b1 = _np.zeros(b2.shape)
            if b2.shape[0] == 0:
                break
        b2, b1 = b1, b2


def match(beamline: _Input,
          beam: _Beam,):
    pass
