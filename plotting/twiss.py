import numpy as np
from georges.plotting.common import palette, filled_plot


def twiss(ax, bl, context, **kwargs):
    """Plot the Twiss beam envelopes from a beamline Twiss computation and a context."""
    bl = bl.line
    ax1 = ax

    bl['XMAXMONO'] = np.sqrt(bl['BETX']*context['EMITX'])
    bl['XMAX'] = np.sqrt(bl['XMAXMONO']**2+(context['DPP']*bl['DX'])**2)
    bl['YMAX'] = np.sqrt(bl['BETY']*context['EMITY'])

    filled_plot(ax1, bl['S'], 0, -1000 * bl['XMAXMONO'], palette['X'], True, alpha=0.8)
    filled_plot(ax1, bl['S'], 0, -1000 * bl['XMAX'], palette['X'], True, alpha=0.4)
    filled_plot(ax1, bl['S'], 0, -2 * 1000 * bl['XMAX'], palette['X'], True, alpha=0.2)
    filled_plot(ax1, bl['S'], 0, 1000 * bl['YMAX'], palette['Y'], True, alpha=0.4)
    filled_plot(ax1, bl['S'], 0, 2 * 1000 * bl['YMAX'], palette['Y'], True, alpha=0.2)


