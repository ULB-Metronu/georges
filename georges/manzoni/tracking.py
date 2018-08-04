import numpy as np
import pandas as pd
from ..beam import Beam
from georges.beamline import Beamline
from . import manzoni
from .common import convert_line
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
            if isinstance(line, _model.ManzoniModel):
                raise TrackException("The line must be a regular Beamline, not a converted Manzoni line.")
            manzoni_line = convert_line(line.line, context)
            manzoni_beam = np.array(beam.distribution)
    else:
        # Access the object's model
        if not isinstance(model, _model.ManzoniModel) and hasattr(model, 'manzoni_model'):
            if hasattr(model, 'model'):
                line = model.model.beamline
            else:
                raise TrackException("The model must also contain a regular Beamline object.")
            model = model.manzoni_model

        elif not isinstance(model, _model.Model) and hasattr(model, 'model'):
            model = model.model

        # Access or convert the beam and beamline to Manzoni's format
        if isinstance(model, _model.ManzoniModel):
            manzoni_line = model.beamline
            manzoni_beam = model.beam
        elif isinstance(model, _model.Model):
            line = model.beamline
            beam = model.beam
            context = model.context
            manzoni_line = convert_line(line.line, context)
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
