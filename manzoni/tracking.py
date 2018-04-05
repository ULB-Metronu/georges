import numpy as np
import pandas as pd
from .manzoni import convert_line
from .observers import ElementByElementObserver
from ..beamline import Beamline
from ..beam import Beam
from . import manzoni


class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def track(line, beam, **Kwargs):
    """Compute the distribution of the beam as it propagates through the beamline.
    """
    # Process arguments
    if line is None or beam is None:
        raise TrackException("Beamline, Beam and MAD-X objects need to be defined.")
    manzoni_line = convert_line(line.line)
    manzoni_beam = np.array(beam.distribution)

    # Run Manzoni
    o = ElementByElementObserver()
    manzoni.track(manzoni_line, manzoni_beam, observer=o)

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

