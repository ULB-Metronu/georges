import os, re, io
import pandas as pd
from .. import beamline
from .. import beam
from .bdsim import BDSim
from .. import physics

class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def read_tracking(file):
    """Read a BDSIM Tracking 'one' file to a dataframe."""
    # TO DO


def track(**kwargs):
    """Compute the distribution of the beam as it propagates through the beamline..
    :param kwargs: parameters are:
        - line: the beamline on which twiss will be run
        - context: the associated context on which MAD-X is run
    """
    # Process arguments
    line = kwargs.get('line', None)
    b = kwargs.get('beam', None)
    context = kwargs.get('context', None)

    if line is None or b is None or context is None:
        raise TrackException("Beamline, Beam, Context and BDsim objects need to be defined.")

    print(context)
    bd = BDSim(beamlines=[line], **kwargs)

    # Create a new beamline to include the results
    l = line.line.copy()

    # Run G4Beamline
    #bd.track(round(g4_beam, 5))
    #errors = bd.run(**kwargs).fatals

    #if kwargs.get("debug", False):
    #    print(g4.raw_input)
    #    print(g4.input)
    #if len(errors) > 0:
    #    print(errors)
    #    # raise TrackException("G4Beamline ended with fatal error.")

    # Add columns which contains datas
    #l['BEAM'] = l.apply(lambda g: read_g4beamline_tracking('Detector' + g.name + '.txt'), axis=1)
    #l.apply(lambda g:
     #       os.remove('Detector' + g.name + '.txt') if os.path.isfile('Detector' + g.name + '.txt') else None,
     #       axis=1)

    #return beamline.Beamline(l)

