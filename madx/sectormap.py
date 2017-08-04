import os
import pandas as pd
from .. import beamline
from .madx import Madx

MADX_SECTORMAP_HEADERS_SKIP_ROWS = 1
MADX_SECTORMAP_DATA_SKIP_ROWS = 8


class SectormapException(Exception):
    """Exception raised for errors in the Sectormap module."""

    def __init__(self, m):
        self.message = m


def read_madx_sectormap(file):
    """Read a MAD-X Sectormap TFS file to a dataframe."""
    headers = pd.read_csv(file, skiprows=MADX_TWISS_HEADERS_SKIP_ROWS, nrows=0, delim_whitespace=True)
    headers.drop(headers.columns[[0, MADX_SECTORMAP_HEADERS_SKIP_ROWS]], inplace=True, axis=1)
    df = pd.read_csv(file,
                     header=None,
                     names=headers,
                     na_filter=False,
                     skiprows=MADX_SECTORMAP_DATA_SKIP_ROWS,
                     delim_whitespace=True
                     )
    df.index.name = 'NAME'
    return df


def sectormap(**kwargs):
    """Compute the transfer matrices of the beamline.
    :param kwargs: parameters are:
        - line: the beamline on which twiss will be run
        - context (optional): the associated context on which MAD-X is run
    """
    # Process arguments
    line = kwargs.get('line', None)
    if line is None:
        raise TwissException("Beamline and MAD-X objects need to be defined.")
    m = Madx(beamlines=line)
    m.beam(line.name)
    m.sectormap(line=kwargs.get('periodic', False),
                start=kwargs.get("start", None),
                places=kwargs.get("places", [])
                )
    errors = m.run(**kwargs).fatals
    if kwargs.get("debug", False):
        print(m.input)
    if len(errors) > 0:
        print(errors)
        raise TwissException("MAD-X ended with fatal error.")
    madx_sectormap = read_madx_sectormap(os.path.join(".", 'sectormap'))
    line_with_sectormap = madx_sectormap.merge(line.line,
                                               left_index=True,
                                               right_index=True,
                                               how='outer',
                                               suffixes=('_SECTORMAP', '')
                                               ).sort_values(by='S')
    return beamline.Beamline(line_with_sectormap)
