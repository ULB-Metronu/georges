import os
import pandas as pd
from .. import beamline
from .madx import Madx

MADX_SECTORMAP_HEADERS_SKIP_ROWS = 6
MADX_SECTORMAP_DATA_SKIP_ROWS = 8


class MatchException(Exception):
    """Exception raised for errors in the Match module."""

    def __init__(self, m):
        self.message = m


def match(**kwargs):
    """Compute the transfer matrices of the beamline.
    :param kwargs: parameters are:
        - line: the beamline on which twiss will be run
        - context: the associated context on which MAD-X is run
        - periodic: ring (True) or beamline (False)
        - start: RANGE at which the twiss/sectormap shoudl start
        - places: where to extract the sector map from
    """
    # Process arguments
    line = kwargs.get('line', None)
    if line is None:
        raise TwissException("Beamline and MAD-X objects need to be defined.")
    m = Madx(beamlines=[line])
    m.beam(line.name)
    m.sectormap(name=line.name,
                line=kwargs.get('periodic', False),
                start=kwargs.get("start", None),
                places=kwargs.get("places", []),
                sectoracc=kwargs.get("sectoracc", False),
                reflect=kwargs.get("reflect", False),
                )
    errors = m.run(**kwargs).fatals
    if kwargs.get("debug", False):
        print(m.input)
    if len(errors) > 0:
        print(errors)
        raise SectormapException("MAD-X ended with fatal error.")
    madx_sectormap = read_madx_sectormap(os.path.join(".", 'sectormap')).rename(columns={'POS': 'S'})
    line_with_sectormap = madx_sectormap.merge(line.line,
                                               left_index=True,
                                               right_index=True,
                                               how='outer',
                                               suffixes=('_SECTORMAP', '')
                                               ).sort_values(by='S')
    return beamline.Beamline(line_with_sectormap)
