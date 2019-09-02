from __future__ import annotations
import numpy as _np
from .input import Input as _Input


def track(beamline: _Input, beam: _np.ndarray):
    b1 = beam  # .copy()
    b2 = _np.zeros(beam.shape)
    for e in beamline.sequence:
        b1, b2 = e.propagate(b1, b2)
        b1, b2 = e.aperture(b1, b2)
        b2 = _np.zeros(b1.shape)


def twiss():
    pass
