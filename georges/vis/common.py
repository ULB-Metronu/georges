import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib.ticker import *

PALETTE = {
    'solarized':  {'base03':  '#002b36',
                   'base02':  '#073642',
                   'base01':  '#586e75',
                   'base00':  '#657b83',
                   'base0':   '#839496',
                   'base1':   '#93a1a1',
                   'base2':   '#eee8d5',
                   'base3':   '#fdf6e3',
                   'yellow':  '#b58900',
                   'orange':  '#cb4b16',
                   'red':     '#dc322f',
                   'magenta': '#d33682',
                   'violet':  '#6c71c4',
                   'blue':    '#268bd2',
                   'cyan':    '#2aa198',
                   'green':   '#859900'
                   }
}

# Define default color palette
palette = PALETTE['solarized']

# Define "logical" colors
palette['quad'] = palette['blue']
palette['bend'] = palette['red']
palette['coll'] = palette['yellow']
palette['X'] = palette['cyan']
palette['Y'] = palette['orange']
palette['X_MADX'] = palette['cyan']
palette['Y_MADX'] = palette['orange']
palette['X_G4BL'] = palette['magenta']
palette['Y_G4BL'] = palette['green']
palette['quadrupole'] = palette['red']
palette['sextupole'] = palette['green']
palette['octupole'] = palette['green']
palette['multipole'] = palette['green']


def style_boxplot(bp, color):
    """Apply fancy styles to a matplotlib boxplot."""
    for box in bp['boxes']:
        box.set(color=color, linewidth=1)
        box.set(facecolor=color, alpha=0.4)
    for whisker in bp['whiskers']:
        whisker.set(color=color, linewidth=1, alpha=0.5)
    for cap in bp['caps']:
        cap.set(color=color, linewidth=1)
    for median in bp['medians']:
        median.set(color=color, linewidth=1)
    for flier in bp['fliers']:
        flier.set(marker='o', markersize=6, color=color, markeredgecolor='none', alpha=0.5)


def beamline_get_ticks_locations(o):
    return list(o['AT_CENTER'].values)


def beamline_get_ticks_labels(o):
    return list(o.index)


def draw_beamline(ax, bl):

    offset = 1.15
    ax2 = ax.twinx()
    ax2.set_yticks([])
    ax2.set_ylim([0, 1])
    ax2.hlines(offset, 0, bl['L'].sum().magnitude, clip_on=False)
    for i, e in bl.query("TYPE=='sbend' or TYPE=='rbend'").iterrows():
        if e['ANGLE'] > 0:
            fc = 'r'
        elif e['ANGLE'] < 0:
            fc = 'b'
        else:
            fc = 'k'
        if e['K1'].magnitude > 0:
            focusing = 1.0
        elif e['K1'].magnitude < 0:
            focusing = -1.0
        else:
            focusing = 0.0
        ax2.add_patch(
                matplotlib.patches.Rectangle(
                    (e['AT_ENTRY'], offset-0.05+focusing*0.02),
                    e['L'].magnitude,
                    .1,
                    hatch='',
                    facecolor='blue',
                    clip_on=False,
                )
            )
    for i, e in bl.query("TYPE=='sextupole' or TYPE=='quadrupole' or TYPE=='multipole'").iterrows():
        fc = 'g'
        ax2.add_patch(
                matplotlib.patches.Rectangle(
                    (e['AT_ENTRY'], offset-0.05),
                    e['L'].magnitude,
                    .1,
                    hatch='',
                    facecolor=palette[e['TYPE']],
                    ec=palette[e['TYPE']],
                    clip_on=False,
                )
            )

    for i, e in bl.query("TYPE=='rcol'").iterrows():
        fc = 'g'
        ax2.add_patch(
                matplotlib.patches.Rectangle(
                    (e['AT_ENTRY'], offset-0.05),
                    e['L'].to('meter').magnitude,
                    .1,
                    hatch='',
                    facecolor='gold',
                    ec='gold',
                    clip_on=False,
                )
            )

    for i, e in bl.query("TYPE=='degrader'").iterrows():
        fc = 'g'
        ax2.add_patch(
                matplotlib.patches.Rectangle(
                    (e['AT_ENTRY'], offset-0.05),
                    e['L'].to('meter').magnitude,
                    .1,
                    hatch='',
                    facecolor='black',
                    ec='black',
                    clip_on=False,
                )
            )
    for i, e in bl.query("CLASS=='CircularCollimator'").iterrows():
        fc = 'g'
        ax2.add_patch(
                matplotlib.patches.Rectangle(
                    (e['AT_ENTRY'], offset-0.05),
                    e['L'].to('meter').magnitude,
                    .1,
                    hatch='',
                    facecolor='grey',
                    ec='grey',
                    clip_on=False,
                )
            )

def prepare(ax, bl, **kwargs):

    bl_short = bl.reset_index()
    bl_short = bl_short[[not a for a in bl_short['NAME'].str.contains("DRIFT")]]
    bl_short = bl_short.set_index("NAME")

    ticks_locations_short = beamline_get_ticks_locations(bl_short)
    ticks_labels_short = beamline_get_ticks_labels(bl_short)
    ticks_locations = beamline_get_ticks_locations(bl)
    ax.tick_params(axis='both', which='major')
    ax.tick_params(axis='x', labelsize=8)
    ax.xaxis.set_major_locator(FixedLocator(ticks_locations_short))

    ax.set_xlim([ticks_locations[0], ticks_locations[-1]])
    ax.get_xaxis().set_tick_params(direction='out')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=-45)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    ax.set_ylim(kwargs.get('ylim', [-60, 60]))
    ax.set_xlabel('s (m)')
    ax.set_ylabel(r'Beam size (mm)')
    ax.grid(False, alpha=0.25)

    if kwargs.get('print_label', True):
        ax2 = ax.twiny()
        ax2.set_xlim([ticks_locations[0], ticks_locations[-1]])
        ax2.get_xaxis().set_tick_params(direction='out')
        ax2.tick_params(axis='both', which='major')
        ax2.tick_params(axis='x', labelsize=8)
        plt.setp(ax2.xaxis.get_majorticklabels(), rotation=-90)
        ax2.xaxis.set_major_formatter(FixedFormatter(ticks_labels_short))
        ax2.xaxis.set_major_locator(FixedLocator(ticks_locations_short))

    if kwargs.get("size_arrows", False):
        ax.set_yticklabels([str(abs(x)) for x in ax.get_yticks()])
        ax.annotate('', xy=(-0.103, 0.97), xytext=(-0.103, 0.75),
                    arrowprops=dict(arrowstyle="->", color='k'), xycoords=ax.transAxes)
        ax.annotate('', xy=(-0.103, 0.25), xycoords='axes fraction', xytext=(-0.103, 0.03),
                    arrowprops=dict(arrowstyle="<-", color='k'))
        ax.text(-0.126, 0.86, "Vertical", fontsize=7, rotation=90, transform=ax.transAxes)
        ax.text(-0.126, 0.22, "Horizontal", fontsize=7, rotation=90, transform=ax.transAxes)

    if kwargs.get("with_beamline", False):
        draw_beamline(ax, bl)


def filled_plot(ax, x, y0, y, c, fill=False, **kwargs):
    ax.plot(x, y, '.', markersize=0,
            markerfacecolor=c, markeredgecolor=c, color=c, **kwargs)
    if fill:
        ax.fill_between(x, y0, y, facecolor=c, linewidth=0.0, edgecolor=c, **kwargs)
