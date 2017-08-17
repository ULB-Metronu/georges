import os
import pandas as pd
from .. import beamline
from .. import beam
from .g4beamline import G4Beamline
import numpy as np

G4BEAMLINE_SKIP_ROWS = 3


class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def read_g4beamline_tracking(file):
    """Read a G4Beamline Tracking 'one' file to a dataframe."""

    column_names = ['X', 'Y', 'S', 'PX', 'PY', 'PZ', 't', 'PDGid', 'EventID', 'TrackID', 'ParentID', 'Weight']
    tmp=np.nan
    if os.path.isfile(file):
        data = pd.read_csv(file, skiprows=G4BEAMLINE_SKIP_ROWS, delimiter=' ', header=None, names=column_names)

        if len(data) == 0:
            return tmp
        data['X'] /= 1000
        data['Y'] /= 1000
        data['S'] /= 1000
        data['P']=np.sqrt(data['PX']**2+data['PY']**2+data['PZ']**2)
        tmp=beam.Beam(data[['X', 'PX', 'Y', 'PY', 'P']])

    return tmp


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
        raise TrackException("Beamline, Beam, context and G4Beamline objects need to be defined.")
    g4 = G4Beamline(beamlines=line, **kwargs)

    # Convert m in mm for G4Beamline
    g4_beam = b.distribution.copy()

    g4_beam['X'] *= 1000
    g4_beam['Y'] *= 1000

    # Create a new beamline to include the results
    l = line.line.copy()

    # Run G4Beamline
    g4.track(g4_beam)
    errors = g4.run(**kwargs).fatals

    if kwargs.get("debug", False):
         print(g4.raw_input)
         print(g4.input)
    if len(errors) > 0:
         print(errors)
         #raise TrackException("G4Beamline ended with fatal error.")

    # Add columns which contains datas
    l['BEAM']=l.apply(lambda g: read_g4beamline_tracking('Detector'+g.name+'.txt'), axis=1)
    l.apply(lambda g:
            os.remove('Detector' + g.name + '.txt') if os.path.isfile('Detector' + g.name + '.txt') else None,
            axis=1)

    return beamline.Beamline(l)
