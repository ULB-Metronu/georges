import os
import pandas as pd
from .. import beamline
from .. import beam
from .madx import Madx

MADX_TRACKING_SKIP_ROWS = 54
PTC_TRACKING_SKIP_ROWS = 9


class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def read_madx_tracking(file):
    """Read a MAD-X Tracking onetable=true file to a dataframe."""
    column_names = ['ID', 'TURN', 'X', 'PX', 'Y', 'PY', 'T', 'PT', 'S', 'E']
    data = pd.read_csv(file, skiprows=MADX_TRACKING_SKIP_ROWS, delim_whitespace=True, names=column_names)
    return data.apply(pd.to_numeric, errors="ignore").dropna()


def read_ptc_tracking(file):
    """Read a PTC Tracking 'one' file to a dataframe."""
    column_names = ['ID', 'TURN', 'X', 'PX', 'Y', 'PY', 'T', 'PT', 'S', 'E']
    data = pd.read_csv(file, skiprows=PTC_TRACKING_SKIP_ROWS, delim_whitespace=True,
                       names=column_names) \
             .apply(pd.to_numeric, errors="ignore").dropna()
    return data[data['TURN'] == 1]


def track(**kwargs):
    """Compute the distribution of the beam as it propagates through the beamline..
    :param kwargs: parameters are:
        - line: the beamline on which twiss will be run
        - context: the associated context on which MAD-X is run
    """
    # Process arguments
    line = kwargs.get('line', None)
    b = kwargs.get('beam', None)
    if line is None or b is None:
        raise TrackException("Beamline, Beam and MAD-X objects need to be defined.")
    m = Madx(beamlines=line)

    # Create a new beamline to include the results
    l = line.line.copy()

    # Run MAD-X
    m.beam(line.name)
    m.track(b.distribution, line, ptc=kwargs.get('ptc', True))
    errors = m.run(**kwargs).fatals
    if kwargs.get("debug", False):
        print(m.raw_input)
        print(m.input)
    if len(errors) > 0:
        print(errors)
        raise TrackException("MAD-X ended with fatal error.")
    if kwargs.get('ptc', True):
        madx_track = read_ptc_tracking(os.path.join(".", 'ptctrackone.tfs'))
    else:
        madx_track = read_madx_tracking(os.path.join(".", 'tracking.outxone')).dropna()
        madx_track['PY'] = pd.to_numeric(madx_track['PY'])
    madx_track['S'] = round(madx_track['S'], 8)
    tmp = madx_track.query('TURN == 1').groupby('S').apply(lambda g: beam.Beam(g[['X', 'PX', 'Y', 'PY', 'PT']]))
    print(tmp)
    l['AT_CENTER_TRUNCATED'] = round(l['AT_CENTER'], 8)
    if 'BEAM' in l:
        l.line.drop('BEAM', inplace=True, axis=1)
    l = l.merge(pd.DataFrame(tmp,
                             columns=['BEAM']),
                left_on='AT_CENTER_TRUNCATED',
                right_index=True,
                how='left').sort_values(by='AT_CENTER')
    l.drop('AT_CENTER_TRUNCATED', axis=1, inplace=True)
    l.sort_values(by='AT_CENTER', inplace=True)
    return beamline.Beamline(l)
