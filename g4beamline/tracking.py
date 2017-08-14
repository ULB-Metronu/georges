import os
import pandas as pd
from .. import beamline
from .. import beam
from .g4beamline import G4Beamline

G4BEAMLINE_SKIP_ROWS = 3


class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def read_g4beamline_tracking(file):
    """Read a G4Beamline Tracking 'one' file to a dataframe."""
    column_names = ['X','Y','S','PX','PY','PZ','t','PDGid','EventID','TrackID','ParentID','Weight']
    data = pd.read_csv(file, skiprows=G4BEAMLINE_SKIP_ROWS,delimiter=' ',header=None,
                       names=column_names)
    return data


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
    g4 = G4Beamline(beamlines=line)

    # Create a new beamline to include the results
    l = line.line.copy()

    # Run G4Beamline
    g4.track(b.distribution)
    errors = g4.run(**kwargs).fatals

    if kwargs.get("debug", False):
         print(g4.raw_input)
         print(g4.input)
    if len(errors) > 0:
         print(errors)
         raise TrackException("G4Beamline ended with fatal error.")

    ## Do the function which reads the file
    # g4_track = read_g4beamline_tracking(os.path.join(".", 'tracking.outxone')).dropna()

    # g4_track['PY'] = pd.to_numeric(madx_track['PY'])
    # g4_track['S'] = round(madx_track['S'], 8)
    # tmp = g4_track.query('TURN == 1').groupby('S').apply(lambda g: beam.Beam(g[['X', 'PX', 'Y', 'PY', 'PT']]))

    # l['AT_CENTER_TRUNCATED'] = round(l['AT_CENTER'], 8)
    # if 'BEAM' in l:
    #     l.line.drop('BEAM', inplace=True, axis=1)
    # l = l.merge(pd.DataFrame(tmp,
    #                          columns=['BEAM']),
    #                          left_on='AT_CENTER_TRUNCATED',
    #                          right_index=True,
    #                          how='left').sort_values(by='AT_CENTER')
    # l.drop('AT_CENTER_TRUNCATED', axis=1, inplace=True)
    return beamline.Beamline(l)
