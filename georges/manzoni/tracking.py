import numpy as np
import pandas as pd
from . import manzoni
from .common import _process_model_argument
from .observers import ElementByElementObserver
from .. import Beamline
from .. import Beam


class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def track(model=None, line=None, beam=None, context={}, **kwargs):
    """Compute the distribution of the beam as it propagates through the beamline.
    """
    # Process arguments
    v = _process_model_argument(model, line, beam, context, TrackException)

    # Run Manzoni
    o = ElementByElementObserver()
    manzoni.track(v['manzoni_line'], v['manzoni_beam'], observer=o, **kwargs)

    # Collect the results
    return Beamline(
        v['georges_line'].line.reset_index().merge(
            pd.DataFrame(
                list(
                    map(
                        lambda x: Beam(
                            pd.DataFrame(x)
                        ),
                        o.data[0, :]
                    )
                ),
                columns=['BEAM']
            ),
            left_index=True,
            right_index=True,
            how='left'
        ).set_index('NAME'))
