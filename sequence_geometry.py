import sys
import inspect
import pandas as pd
import numpy as np


def compute_derived_data(row):
    """Compute derived data for a beamline element."""
    # Add a bunch of mandatory columns
    for c in ['E1', 'E2', 'FINT', 'ANGLE', 'HGAP', 'THICK', 'TILT']:
        if pd.isnull(row.get(c)):
            row[c] = np.nan

    # Corner case
    if pd.isnull(row.get('CLASS')):
        row['CLASS'] = row.get('TYPE', 'MARKER')

    # Should be changed to infer the list of the module functions
    all_converters = [e[0].upper() for e in inspect.getmembers(sys.modules[__name__], inspect.isfunction)
                      if e[0] is not "compute_derived_data"]
    for d in row.index | all_converters:
        if pd.isnull(row.get(d)) and not pd.isnull(row.get(d + '_ELEMENT')):
            row[d] = row[d + '_ELEMENT']
        if pd.isnull(row.get(d)):
            row[d] = getattr(sys.modules[__name__], d.lower(), lambda r: np.nan)(row)

    return row


def at_entry(r):
    """Try to compute the element's entry 's' position from other data."""
    if not pd.isnull(r.get('AT_CENTER')):
        return r['AT_CENTER'] - r['ORBIT_LENGTH']/2.0
    elif not pd.isnull(r.get('AT_EXIT')):
        return r['AT_EXIT'] - r['ORBIT_LENGTH']
    else:
        return np.nan


def at_center(r):
    """Try to compute the element's center 's' position from other data."""
    if pd.isnull(r.get('ORBIT_LENGTH')):
        return np.nan
    if not pd.isnull(r.get('AT_ENTRY')):
        return r['AT_ENTRY'] + r['ORBIT_LENGTH']/2.0
    elif not pd.isnull(r.get('AT_EXIT')):
        return r['AT_EXIT'] - r['ORBIT_LENGTH']/2.0
    else:
        return np.nan


def at_exit(r):
    """Try to compute the element's entry 's' position from other data."""
    if pd.isnull(r.get('ORBIT_LENGTH')):
        return np.nan
    if not pd.isnull(r.get('AT_ENTRY')):
        return r['AT_ENTRY'] + r['ORBIT_LENGTH']
    elif not pd.isnull(r.get('AT_CENTER')):
        return r['AT_CENTER'] + r['ORBIT_LENGTH']/2.0
    else:
        return np.nan


def orbit_length(r):
    """Try to compute the element's orbit length from other data."""
    if pd.isnull(r.get('LENGTH')):
        return 0.0
    if pd.isnull(r.get('ANGLE')):
        return r['LENGTH']
    return r['ANGLE']*np.pi/180.0 * r['LENGTH'] / (2.0 * np.sin((r['ANGLE']*np.pi/180.0) / 2.0))


def truc(r):
    return 1.0