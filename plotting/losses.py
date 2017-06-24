from matplotlib.ticker import *
from .common import beamline_get_ticks_locations
from matplotlib.ticker import *
from .common import *


def losses(ax, transmission, bl):
    ticks_locations = beamline_get_ticks_locations(bl.line)
    ax2 = ax.twinx()

    ax2.get_xaxis().set_tick_params(direction='out')
    ax2.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('right')
    ax2.set_xlim([ticks_locations[0], ticks_locations[-1]])
    ax2.tick_params(axis='x', labelsize=6)
    ax2.xaxis.set_major_formatter(FixedFormatter([]))
    ax2.xaxis.set_major_locator(FixedLocator(ticks_locations))
    ax2.yaxis.set_major_locator(MultipleLocator(25))
    ax2.set_ylabel('T ($\%$)')
    ax2.yaxis.label.set_color(palette['green'])
    ax2.set_ylim([0, 100])
    ax2.grid(True)
    ax2.plot(transmission, '^-', color=palette['green'])

    ax.set_xlim([ticks_locations[0], ticks_locations[-1]])
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.set_ylabel('Losses ($\%$)')
    ax.yaxis.label.set_color(palette['magenta'])
    ax.bar(transmission.index-0.125, -transmission.diff(), 0.125,alpha=0.7,
            edgecolor=palette['magenta'],
            color=palette['magenta'],
            #yerr=transmission.apply(compute_losses_error),
            error_kw=dict(ecolor=palette['base02'], lw=1, capsize=2, capthick=1))
    ax.set_ylim([0, transmission.diff().abs().max()+5.0])
