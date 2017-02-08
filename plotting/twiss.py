import numpy as np
from georges.plotting.common import palette, filled_plot


def twiss(ax, twiss, context, **kwargs):
    ax1 = ax

    twiss['XMAXMONO'] = np.sqrt(twiss['BETX']*context['EMITX'])
    twiss['XMAX'] = np.sqrt(twiss['XMAXMONO']**2+(context['DPP']*twiss['DX'])**2)
    twiss['YMAX'] = np.sqrt(twiss['BETY']*context['EMITY'])

    filled_plot(ax1, twiss['S'], 0, -1000 * twiss['XMAXMONO'], palette['X'], True, alpha=0.8)
    filled_plot(ax1, twiss['S'], 0, -1000 * twiss['XMAX'], palette['X'], True, alpha=0.4)
    filled_plot(ax1, twiss['S'], 0, -2 * 1000 * twiss['XMAX'], palette['X'], True, alpha=0.2)
    filled_plot(ax1, twiss['S'], 0, 1000 * twiss['YMAX'], palette['Y'], True, alpha=0.4)
    filled_plot(ax1, twiss['S'], 0, 2 * 1000 * twiss['YMAX'], palette['Y'], True, alpha=0.2)


