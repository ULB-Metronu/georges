import numpy as np
from .common import palette, filled_plot


def twiss(ax, bl, context, **kwargs):
    """Plot the Twiss beam envelopes from a beamline Twiss computation and a context."""
    bl = bl.line

    bl['XMAXMONO'] = np.sqrt(bl['BETX']*context['EMITX'])
    bl['XMAX'] = np.sqrt(bl['XMAXMONO']**2+(context['DPP']*bl['DX'])**2)
    bl['YMAX'] = np.sqrt(bl['BETY']*context['EMITY'])

    filled_plot(ax, bl['S'], 0, -1000 * bl['XMAXMONO'], palette['X'], True, alpha=0.8)
    filled_plot(ax, bl['S'], 0, -1000 * bl['XMAX'], palette['X'], True, alpha=0.4)
    filled_plot(ax, bl['S'], 0, -2 * 1000 * bl['XMAX'], palette['X'], True, alpha=0.2)
    filled_plot(ax, bl['S'], 0, 1000 * bl['YMAX'], palette['Y'], True, alpha=0.4)
    filled_plot(ax, bl['S'], 0, 2 * 1000 * bl['YMAX'], palette['Y'], True, alpha=0.2)


def beta(ax, bl, **kwargs):
    """Plot the Twiss beta functions."""
    bl = bl.line

    plt.plot(ax, bl['S'], bl['BETX'])
    plt.plot(ax, bl['S'], bl['BETY'])


def alpha(ax, bl, context, **kwargs):
    """Plot the Twiss alpha functions."""
    bl = bl.line

    plt.plot(ax, bl['S'], bl['ALFX'])
    plt.plot(ax, bl['S'], bl['ALFY'])


def dispersion(ax, bl, **kwargs):
    """Plot the dispersion functions."""
    bl = bl.line

    plt.plot(ax, bl['S'], bl['DISPX'])
    plt.plot(ax, bl['S'], bl['DISPY'])