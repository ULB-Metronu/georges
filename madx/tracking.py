import os
import pandas as pd
import georges.beamline as beamline
import georges.beam as beam


class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def read_madx_tracking(file):
    """Read a MAD-X Tracking onetable=true file to a dataframe."""
    column_names = ['ID','TURN','X','PX','Y','PY','T','PT','S','E']
    data = pd.read_csv(file, skiprows=54, delim_whitespace=True, names=column_names)
    return data.apply(pd.to_numeric, errors="ignore").dropna()


def read_ptc_tracking(file):
    """Read a PTC Tracking 'one' file to a dataframe."""
    column_names = ['ID', 'TURN', 'X', 'PX', 'Y', 'PY', 'T', 'PT', 'S', 'E']
    data = pd.read_csv(file, skiprows=9, delim_whitespace=True,
                       names=column_names) \
              .apply(pd.to_numeric, errors="ignore").dropna()
    return data[data['TURN'] == 1]


def track(**kwargs):
    """Compute the distribution of the beam as it propagates through the beamline."""
    # Process arguments
    line = kwargs.get('line', None)
    m = kwargs.get('madx', None)
    b = kwargs.get('beam', None)
    if line is None or m is None or b is None:
        raise TrackException("Beamline, Beam and MAD-X objects need to be defined.")

    # Create a new beamline to include the results
    l = line.line.copy()

    # Attach the new beamline to MAD-X if needed
    if line not in m.beamlines:
        m.attach(line)
    m.beam(line.name)
    m.track(b.distribution, line, ptc=kwargs.get('ptc', True))
    errors = m.run(m.context).fatals
    if len(errors) > 0:
        m.print_input()
        print(errors)
        raise TrackException("MAD-X ended with fatal error.")
    if kwargs.get('ptc', True):
        madx_track = read_ptc_tracking(os.path.join(m.path, 'ptctrackone.tfs'))
    else:
        madx_track = read_madx_tracking(os.path.join(m.path, 'tracking.outxone')).dropna()
        madx_track['PY'] = pd.to_numeric(madx_track['PY'])
    madx_track['S'] = round(madx_track['S'], 8)
    tmp = madx_track.query('TURN == 1').groupby('S').apply(lambda g: beam.Beam(g[['X', 'PX', 'Y', 'PY', 'PT']]))
    l['AT_CENTER_TRUNCATED'] = round(l['AT_CENTER'], 8)
    if 'BEAM' in l:
        l.line.drop('BEAM', inplace=True, axis=1)
    print(l)
    l = l.merge(pd.DataFrame(tmp, columns=['BEAM']),
                             left_on='AT_CENTER_TRUNCATED',
                             right_index=True,
                             how='left'
                          ).sort_values(by='AT_CENTER')
    l.drop('AT_CENTER_TRUNCATED', axis=1, inplace=True)
    l.sort_values(by='S', inplace=True)
    return beamline.Beamline(l)
