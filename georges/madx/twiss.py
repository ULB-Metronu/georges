import os
import re
import pandas as pd
from .. import beamline
from .madx import Madx

MADX_TWISS_HEADERS_SKIP_ROWS = 45
MADX_TWISS_DATA_SKIP_ROWS = 47
PTC_TWISS_HEADERS_SKIP_ROWS = 88
PTC_TWISS_DATA_SKIP_ROWS = 90


class TwissException(Exception):
    """Exception raised for errors in the Twiss module."""

    def __init__(self, m):
        self.message = m


def read_madx_sectormap(file):
    """Read a MAD-X Sectormap TFS file to a dataframe."""
    headers = pd.read_csv(file, skiprows=MADX_TWISS_HEADERS_SKIP_ROWS, nrows=0, delim_whitespace=True)
    headers.drop(headers.columns[[0, 1]], inplace=True, axis=1)
    df = pd.read_csv(file,
                     header=None,
                     names=headers.columns.values,
                     na_filter=False,
                     skiprows=MADX_TWISS_DATA_SKIP_ROWS,
                     delim_whitespace=True
                     )
    df.index.name = 'NAME'
    return df


def read_madx_twiss(file):
    """Read a MAD-X Twiss TFS file to a dataframe."""
    headers = pd.read_csv(file, skiprows=MADX_TWISS_HEADERS_SKIP_ROWS, nrows=0, delim_whitespace=True)
    headers.drop(headers.columns[[0, 1]], inplace=True, axis=1)
    df = pd.read_csv(file,
                     header=None,
                     names=headers.columns.values,
                     na_filter=False,
                     skiprows=MADX_TWISS_DATA_SKIP_ROWS,
                     delim_whitespace=True
                     )
    df.index.name = 'NAME'
    return df


def read_ptc_twiss(file):
    """Read a MAD-X PTC Twiss TFS file to a dataframe."""
    headers = pd.read_csv(file, skiprows=PTC_TWISS_HEADERS_SKIP_ROWS, nrows=0, delim_whitespace=True)
    headers.drop(headers.columns[[0, 1]], inplace=True, axis=1)
    df = pd.read_csv(file,
                     header=None,
                     names=headers.columns.values,
                     na_filter=False,
                     skiprows=PTC_TWISS_DATA_SKIP_ROWS,
                     delim_whitespace=True
                     )
    df.index.name = 'NAME'
    return df


def read_twiss_summary(file):
    """Read a MAD-X PTCS Twiss TFS file summary information (summary table) to a dictionary."""
    regex_header = re.compile("^@\s+(\S*)\s+\S+\s+(.*)$")
    summary = {}
    for line in open("/Users/chernals/Google Drive/Stage IBA/2018/Notebooks/ptc_twiss.outx"):
        matches = re.findall(regex_header, line)
        if len(matches) >= 1:
            matches = matches[0]
        else:
            break
        try:
            summary[matches[0]] = float(matches[1])
        except ValueError:
            summary[matches[0]] = matches[1]
    return summary


def twiss(with_summary=False, **kwargs):
    """Compute the Twiss parameters of the beamline.
    :param with_summary: activate the reading of the MAD-X/PTC summary table.
    :param kwargs: parameters are:
        - line: the beamline on which twiss will be run
        - context: the associated context on which MAD-X is run
    """
    # Process arguments
    line = kwargs.get('line', None)
    if line is None:
        raise TwissException("Beamline and MAD-X objects need to be defined.")
    m = Madx(beamlines=[line], ptc_use_knl_only=kwargs.get('ptc_use_knl_only', False), context=kwargs.get('context'))
    m.beam(line.name)
    m.twiss(**kwargs)
    errors = m.run(**kwargs).fatals
    if kwargs.get("debug", False):
        print(m.input)
    if len(errors) > 0:
        print(errors)
        raise TwissException("MAD-X ended with fatal error.")
    if kwargs.get('ptc', False):
        madx_twiss = read_ptc_twiss(os.path.join(".", 'ptc_twiss.outx'))
    else:
        madx_twiss = read_madx_twiss(os.path.join(".", 'twiss.outx'))
    line_with_twiss = madx_twiss.merge(line.line,
                                       left_index=True,
                                       right_index=True,
                                       how='outer',
                                       suffixes=('_TWISS', '')
                                       ).sort_values(by='S')
    if with_summary:
        madx_twiss_summary = read_twiss_summary(os.path.join(".", 'ptc_twiss.outx'))
        return {
            'line': beamline.Beamline(line_with_twiss),
            'summary': madx_twiss_summary,
        }
    else:
        return beamline.Beamline(line_with_twiss)
