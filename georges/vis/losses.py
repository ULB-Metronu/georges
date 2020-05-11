from .common import beamline_get_ticks_locations
from .common import palette as common_palette
from .common import prepare
from matplotlib.ticker import *
import pandas as pd


def losses(ax, bl,beam_o_df, **kwargs):
    """Plot the losses from a beamline tracking computation and a context."""
    palette = kwargs.get("palette", common_palette)
    bl['n_particles'] = list(map(len, beam_o_df['BEAM_OUT']))
    init = bl.query("TYPE == TYPE").drop_duplicates(subset='AT_CENTER', keep='first').iloc[0]['n_particles']
    transmission = bl.query("TYPE == TYPE").drop_duplicates(subset='AT_CENTER', keep='first').apply(
        lambda r: pd.Series({
            'S': r['AT_EXIT'],
            'T': r['n_particles'] / init
        }), axis=1)

    ax2 = ax.twinx()
    prepare(ax2,bl)
    prepare(ax, bl, with_beamline=kwargs.get("with_beamline", False))

    ticks_locations = beamline_get_ticks_locations(bl)
    #ax2.get_xaxis().set_tick_params(direction='out')
    #ax2.yaxis.set_ticks_position('left')
    #ax2.xaxis.set_ticks_position('bottom')
    #ax2.yaxis.set_ticks_position('right')
    #ax2.set_xlim([ticks_locations[0], ticks_locations[-1]])
    #ax2.tick_params(axis='x', labelsize=6)
    #ax2.xaxis.set_major_formatter(FixedFormatter([]))
    #ax2.xaxis.set_major_locator(FixedLocator(ticks_locations))
    ax2.set_ylabel('T ($\%$)')
    ax2.yaxis.label.set_color(palette['green'])
    ax2.grid(True)
    if kwargs.get('log', False):
        ax2.semilogy(transmission['S'], 100*transmission['T'], 's-', color=palette['green'])
    else:
        ax2.yaxis.set_major_locator(MultipleLocator(10))
        ax2.set_ylim([0, 100])
        ax2.plot(transmission['S'], 100*transmission['T'], 's-', color=palette['green'])
    ax.set_xlim([ticks_locations[0], ticks_locations[-1]])
    #ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.set_ylabel('Losses ($\%$)')
    ax.yaxis.label.set_color(palette['magenta'])
    ax.bar(transmission['S'] - 0.125, -100 * transmission['T'].diff(), 0.125, alpha=0.7,
           edgecolor=palette['magenta'],
           color=palette['magenta'],
           error_kw=dict(ecolor=palette['base02'], lw=1, capsize=2, capthick=1))
    ax.set_ylim([0, 100 * transmission['T'].diff().abs().max() + 5.0])

    if kwargs.get('with_current'):
        ax2.plot(bl['AT_EXIT'], bl['CURRENT'], 'k-')

    return transmission
