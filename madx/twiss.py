import os
import pandas as pd
import georges.beamline as beamline


class TwissException(Exception):
    """Exception raised for errors in the Twiss module."""

    def __init__(self, m):
        self.message = m


def read_madx_twiss(file):
    """Read a MAD-X Twiss TFS file to a dataframe."""
    headers = pd.read_csv(file, skiprows=45, nrows=0, delim_whitespace=True)
    headers.drop(headers.columns[[0,1]], inplace=True, axis=1)
    df = pd.read_csv(file, header=None, names=headers, na_filter=False, skiprows=47, delim_whitespace=True)
    df.index.name = 'NAME'
    return df


def read_ptc_twiss(file):
    """Read a MAD-X PTC Twiss TFS file to a dataframe."""
    headers = pd.read_csv(file, skiprows=88, nrows=0, delim_whitespace=True)
    headers.drop(headers.columns[[0,1]], inplace=True, axis=1)
    df = pd.read_csv(file, header=None, names=headers, na_filter=False, skiprows=90, delim_whitespace=True)
    df.index.name = 'NAME'
    return df


def twiss(**kwargs):
    """Compute the Twiss parameters of the beamline."""
    # Process arguments
    line = kwargs.get('line', None)
    m = kwargs.get('madx', None)
    if line is None or m is None:
        raise TwissException("Beamline and MAD-X objects need to be defined.")
    # Attach the new beamline to MAD-X if needed
    if line not in m.beamlines:
        m.attach(line)
    m.beam(line.name)

    m.twiss(ptc=kwargs.get('ptc', False), centre=True)
    errors = m.run(m.context).fatals
    if len(errors) > 0:
        m.print_input()
        print(errors)
        raise TwissException("MAD-X ended with fatal error.")
    if kwargs.get('ptc', False):
        madx_twiss = read_ptc_twiss(os.path.join(m.path, 'ptc_twiss.outx'))
    else:
        madx_twiss = read_madx_twiss(os.path.join(m.path, 'twiss.outx'))
    line_with_twiss = madx_twiss.merge(line.line,
                                       left_index=True,
                                       right_index=True,
                                       how='outer',
                                       suffixes=('_TWISS', '')
                                       ).sort_values(by='S')
    return beamline.Beamline(line_with_twiss)
