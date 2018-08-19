import numpy as np
import pandas as pd
from ..beam import Beam
from georges.beamline import Beamline
from . import manzoni
from .common import convert_line
from .observers import ElementByElementObserver
from .. import model as _model


class SigmaException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def sigma(model=None, line=None, beam=None, context={}, **kwargs):
    """Compute the distribution of the beam as it propagates through the beamline.
    """
    # Process arguments


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
