import numpy as np
import pandas as pd
from ..beam import Beam
from georges.beamline import Beamline
from . import manzoni
from .manzoni import convert_line
from .observers import ElementByElementObserver


class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def track(model=None, line=None, beam=None, **kwargs):
    """Compute the distribution of the beam as it propagates through the beamline.
    """
    # Process arguments
    if model is None:
        if line is None or beam is None:
            raise TrackException("Beamline and Beam objects need to be defined.")
    else:
        line = model.beamline
        beam = model.beam
    manzoni_line = convert_line(line.line, kwargs.get('context', {}))
    manzoni_beam = np.array(beam.distribution)

    # Run Manzoni
    o = ElementByElementObserver()
    manzoni.track(manzoni_line, manzoni_beam, observer=o, **kwargs)

    # Collect the results
    return Beamline(
        line.line.reset_index().merge(
            pd.DataFrame(
                list(
                    map(
                        lambda x: Beam(
                            pd.DataFrame(x)
                        ),
                        o.data[0, :, 0]
                    )
                ),
                columns=['BEAM']
            ),
            left_index=True,
            right_index=True,
            how='left'
        ).set_index('NAME'))

