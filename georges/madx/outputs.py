from typing import List
import os
import pandas as _pd
from .. import ureg as _ureg

MADX_TWISS_TABLE_HEADER_ROWS: int = 47
"""MAD-X Twiss table header rows (lines to be skipped when reading the table's content."""

MADX_TWISS_HEADERS: List[str] = [
    'NAME',
    'KEYWORD',
    'S',
    'BETX', 'ALFX', 'MUX',
    'BETY', 'ALFY', 'MUY',
    'X', 'PX', 'Y', 'PY', 'T', 'PT',
    'DX', 'DPX', 'DY', 'DPY',
    'WX',
    'PHIX',
    'DMUX',
    'WY',
    'PHIY',
    'DMUY',
    'DDX', 'DDPX',
    'DDY', 'DDPY',
    'R11', 'R12', 'R21', 'R22',
    'ENERGY',
    'L',
    'ANGLE',
    'K0L', 'K0SL',
    'K1L', 'K1SL',
    'K2L', 'K2SL',
    'K3L', 'K3SL',
    'K4L', 'K4SL',
    'K5L', 'K5SL',
    'K6L', 'K6SL',
    'K7L', 'K7SL',
    'K8L', 'K8SL',
    'K9L', 'K9SL',
    'K10L', 'K10SL',
    'K11L', 'K11SL',
    'K12L', 'K12SL',
    'K13L', 'K13SL',
    'K14L', 'K14SL',
    'K15L', 'K15SL',
    'K16L', 'K16SL',
    'K17L', 'K17SL',
    'K18L', 'K18SL',
    'K19L', 'K19SL',
    'K20L', 'K20SL',
    'KSI',
    'HKICK', 'VKICK',
    'TILT',
    'E1', 'E2',
    'H1', 'H2',
    'HGAP',
    'FINT',
    'FINTX',
    'VOLT',
    'LAG',
    'FREQ',
    'HARMON',
    'SLOT_ID',
    'ASSEMBLY_ID',
    'MECH_SEP',
    'V_POS',
    'BBCHARGE', 'XMA', 'YMA', 'SIGX', 'SIGY',
    'LRAD',
    'PARENT',
    'COMMENTS',
    'RE11', 'RE12', 'RE13', 'RE14', 'RE15', 'RE16',
    'RE21', 'RE22', 'RE23', 'RE24', 'RE25', 'RE26',
    'RE31', 'RE32', 'RE33', 'RE34', 'RE35', 'RE36',
    'RE41', 'RE42', 'RE43', 'RE44', 'RE45', 'RE46',
    'RE51', 'RE52', 'RE53', 'RE54', 'RE55', 'RE56',
    'RE61', 'RE62', 'RE63', 'RE64', 'RE65', 'RE66',
    'KMAX',
    'KMIN',
    'CALIB',
    'POLARITY',
    'ALFA',
    'BETA11', 'BETA12', 'BETA13',
    'BETA21', 'BETA22', 'BETA23',
    'BETA31', 'BETA32', 'BETA33',
    'ALFA11', 'ALFA12', 'ALFA13',
    'ALFA21', 'ALFA22', 'ALFA23',
    'ALFA31', 'ALFA32', 'ALFA33',
    'GAMA11', 'GAMA12', 'GAMA13',
    'GAMA21', 'GAMA22', 'GAMA23',
    'GAMA31', 'GAMA32', 'GAMA33',
    'BETA11P', 'BETA12P', 'BETA13P',
    'BETA21P', 'BETA22P', 'BETA23P',
    'BETA31P', 'BETA32P', 'BETA33P',
    'ALFA11P', 'ALFA12P', 'ALFA13P',
    'ALFA21P', 'ALFA22P', 'ALFA23P',
    'ALFA31P', 'ALFA32P', 'ALFA33P',
    'GAMA11P', 'GAMA12P', 'GAMA13P',
    'GAMA21P', 'GAMA22P', 'GAMA23P',
    'GAMA31P', 'GAMA32P', 'GAMA33P',
    'DISP1', 'DISP2', 'DISP3', 'DISP4',
    'DISP1P', 'DISP2P', 'DISP3P', 'DISP4P',
    'DISP1P2', 'DISP2P2', 'DISP3P2', 'DISP4P2',
    'DISP1P3', 'DISP2P3', 'DISP3P3', 'DISP4P3',
    'MU1', 'MU2', 'MU3',
    'SIG11', 'SIG12', 'SIG13', 'SIG14', 'SIG15', 'SIG16',
    'SIG21', 'SIG22', 'SIG23', 'SIG24', 'SIG25', 'SIG26',
    'SIG31', 'SIG32', 'SIG33', 'SIG34', 'SIG35', 'SIG36',
    'SIG41', 'SIG42', 'SIG43', 'SIG44', 'SIG45', 'SIG46',
    'SIG51', 'SIG52', 'SIG53', 'SIG54', 'SIG55', 'SIG56',
    'SIG61', 'SIG62', 'SIG63', 'SIG64', 'SIG65', 'SIG66',
    'N1',
]
"""MAD-X Twiss headers (by default, when all columns are selected)."""


def load_madx_twiss_headers(filename: str = 'twiss.outx', path: str = '.') -> _pd.Series:
    """

    Args:
        filename: name of the Twiss table file
        path: path to the Twiss table file

    Returns:

    """
    return _pd.read_csv(os.path.join(path, filename),
                        sep=r'\s+',
                        usecols=['KEY', 'VALUE'],
                        squeeze=True,
                        index_col=0,
                        names=['@', 'KEY', '_', 'VALUE'],
                        converters={'PC': float},
                        )[0:46]


def load_madx_twiss_table(filename: str = 'twiss.outx',
                          path: str = '.',
                          columns: List = None,
                          with_units: bool = True,
                          ) -> _pd.DataFrame:
    """

    Args:
        filename: name of the Twiss table file
        path: path to the Twiss table file
        columns: the list of columns in the Twiss file
        with_units:

    Returns:
        A DataFrame representing the Twiss table.
    """
    columns = columns or MADX_TWISS_HEADERS
    _: _pd.DataFrame = _pd \
        .read_csv(os.path.join(path, filename),
                  skiprows=MADX_TWISS_TABLE_HEADER_ROWS,
                  sep=r'\s+',
                  index_col=False,
                  names=columns,
                  ) \
        .drop(0)
    for c in _.columns:
        try:
            _[c] = _[c].apply(float)
        except ValueError:
            pass
    if with_units:
        _['L'] = _['L'].apply(lambda e: e * _ureg.m)
        _['E1'] = _['E1'].apply(lambda e: e * _ureg.radian)
        _['E2'] = _['E2'].apply(lambda e: e * _ureg.radian)
        _['ANGLE'] = _['ANGLE'].apply(lambda e: e * _ureg.radian)
        _['K1L'] = _['K1L'].apply(lambda e: e / _ureg.m)
        _['TILT'] = _['TILT'].apply(lambda e: e * _ureg.radian)
    return _.set_index('NAME')


def get_twiss_values(table: _pd.DataFrame, location: int = 0) -> _pd.Series:
    """Extract the initial Twiss parameters from a Twiss table

    Args:
        table: a MAD-X twiss table read as a DataFrame
        location: the location at which the parameters need to be extracted

    Returns:
        A Pandas Series containing the extracted Twiss parameters.
    """
    return _pd.Series({
        'MU1': 0,
        'MU2': 0,
        'BETA11': table.iloc[location]['BETX'],
        'BETA22': table.iloc[location]['BETY'],
        'ALPHA11': table.iloc[location]['ALFX'],
        'ALPHA22': table.iloc[location]['ALFY'],
        'GAMMA11': (1 + table.iloc[location]['ALFX'] ** 2) / table.iloc[location]['BETX'],
        'GAMMA22': (1 + table.iloc[location]['ALFY'] ** 2) / table.iloc[location]['BETY'],
        'DY': table.iloc[location]['DX'],
        'DYP': table.iloc[location]['DPX'],
        'DZ': table.iloc[location]['DY'],
        'DZP': table.iloc[location]['DPY'],
    })