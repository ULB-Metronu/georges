import numpy as np
import pandas as pd
from ..beam import Beam
from georges.beamline import Beamline
from . import manzoni
from .common import convert_line
from .constants import *
from .matrices import matrices4
from .. import model as _model


class TwissException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


class TwissMap:
    def __init__(self, m):
        if m.shape != (4, 4):
            raise TwissException("Map dimensions must be 4x4.")
        self._m = m

    @property
    def mh(self):
        return self._m[0:2, 0:2]

    @property
    def mv(self):
        return self._m[2:4, 2:4]

    @property
    def stable(self):
        return self.stable_h & self.stable_v

    @property
    def stable_h(self):
        return self._stability(self.mh)

    @property
    def stable_v(self):
        return self._stability(self.mv)

    @staticmethod
    def _stability(m):
        return np.abs(np.trace(m)) < 2

    @property
    def tune_h(self):
        return self._tune(self.mh)

    @property
    def tune_v(self):
        return self._tune(self.mv)

    @staticmethod
    def _tune(m):
        if TwissMap._stability(m):
            return np.arccos(0.5 * np.trace(m)) / (2 * np.pi)
        else:
            return np.arccosh(0.5 * np.trace(m)) / (2 * np.pi)

    @property
    def beta_h(self):
        return self._beta(self.mh)

    @property
    def beta_v(self):
        return self._beta(self.mv)

    @staticmethod
    def _beta(m):
        return m[0, 1] * np.sin(2*np.pi*TwissMap._tune(m))

    @property
    def alpha_h(self):
        return self._alpha(self.mh)

    @property
    def alpha_v(self):
        return self._alpha(self.mv)

    @staticmethod
    def _alpha(m):
        return (m[0, 0] - m[1, 1]) / (2 * np.sin(2*np.pi*TwissMap._tune(m)))

    @property
    def twiss_h(self):
        return {
            'beta': self.beta_h,
            'alpha': self.alpha_h,
        }

    @property
    def twiss_v(self):
        return {
            'beta': self.beta_v,
            'alpha': self.alpha_v,
        }


def compute_map(line, **kwargs):
    """
    Sigma-matrix tracking through a beamline.
    Code optimized for performance.
    :param line: beamline description in Manzoni format
    :param kwargs: optional parameters
    :return: Observer.track_end() return value
    """
    m = np.identity(4)

    # Main loop
    for i in range(0, line.shape[0]):
        if line[i, INDEX_CLASS_CODE] in CLASS_CODE_MATRIX:
            m = np.matmul(matrices4[int(line[i, INDEX_CLASS_CODE])](line[i]), m)
        else:
            # Default to a drift
            m = np.matmul(matrices4[CLASS_CODES['DRIFT']](line[i]), m)
    return m


def compute_twiss(m, alpha, beta, mu):
    """
    Compute Twiss parameters as they are propagated by a transfer matrix.
    :param m: transfer matrix
    :param alpha: alpha parameter at the end of the element
    :param beta: beta parameter at the end of the element
    :param mu: phase advance at the end of the element
    :return:
    """
    return {
        'beta': (1/beta) * ((m[0, 0] * beta - m[0, 1] * alpha)**2 + m[0, 1]**2),
        'alpha': -(1/beta) * (
            (m[0, 0] * beta - m[0, 1] * alpha) * (m[1, 0] * beta - m[1, 1] * alpha) + m[0, 1] * m[1, 0]
        ),
        'mu': (2 * np.pi * mu + np.arctan(m[0, 1] / (m[0, 0] * beta - m[0, 1] * alpha))) / (2 * np.pi)
    }


def twiss(model=None, line=None, context={}, periodic=True, **kwargs):
    """
    """
    # Process arguments
    l = convert_line(line.line)
    line = line.line.copy()
    line['BETA11'] = np.nan
    line['BETA22'] = np.nan
    line['ALPHA11'] = np.nan
    line['ALPHA22'] = np.nan
    line['MU1'] = np.nan
    line['MU2'] = np.nan

    # Do the one-turn map computation
    twiss_map = TwissMap(compute_map(l))
    # line.line.iloc[0]['BETA11'] = twiss_map.beta_h
    # line.line.iloc[0]['BETA22'] = twiss_map.beta_v
    # line.line.iloc[0]['ALPHA11'] = twiss_map.alpha_h
    # line.line.iloc[0]['ALPHA22'] = twiss_map.alpha_v
    # line.line.iloc[0]['MU1'] = 0
    # line.line.iloc[0]['MU2'] = 0

    # Do the twiss computation
    twiss_h = twiss_map.twiss_h
    twiss_h['mu'] = 0
    twiss_v = twiss_map.twiss_v
    twiss_v['mu'] = 0
    for i in range(0, l.shape[0]):
        if l[i, INDEX_CLASS_CODE] in CLASS_CODE_MATRIX:
            m = matrices4[int(l[i, INDEX_CLASS_CODE])](l[i])
        else:
            # Default to a drift
            m = matrices4[CLASS_CODES['DRIFT']](l[i])
        twiss_h = compute_twiss(m[0:2, 0:2], **twiss_h)
        twiss_v = compute_twiss(m[2:4, 2:4], **twiss_v)
        line.iat[i, line.columns.get_loc('BETA11')] = twiss_h['beta']
        line.iat[i, line.columns.get_loc('BETA22')] = twiss_v['beta']
        line.iat[i, line.columns.get_loc('ALPHA11')] = twiss_h['alpha']
        line.iat[i, line.columns.get_loc('ALPHA22')] = twiss_v['alpha']
        line.iat[i, line.columns.get_loc('MU1')] = twiss_h['mu']
        line.iat[i, line.columns.get_loc('MU2')] = twiss_v['mu']

    return {
        'map': twiss_map,
        'line': line,
    }
