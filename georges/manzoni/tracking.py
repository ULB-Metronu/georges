import numpy as np
import pandas as pd
from ..beam import Beam
from georges.beamline import Beamline
from . import manzoni
from .manzoni import convert_line
from .observers import ElementByElementObserver
from .. import model as _model


class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def track(model=None, line=None, beam=None, context={}, **kwargs):
    """Compute the distribution of the beam as it propagates through the beamline.
    """
    # Process arguments
    if model is None:
        if line is None or beam is None:
            raise TrackException("Beamline and Beam objects need to be defined.")
        else:
            manzoni_line = convert_line(line.line, context)
            manzoni_beam = np.array(beam.distribution)
    else:
        if not isinstance(model, _model.Model) and hasattr(model, 'model'):
            model = model.model
            line = model.beamline
            beam = model.beam
            context = model.context
            manzoni_line = convert_line(line.line, context)
            manzoni_beam = np.array(beam.distribution)
        elif isinstance(model, model.ManzoniModel):
            manzoni_line = model.beamline
            manzoni_beam = model.beam

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

