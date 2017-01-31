import re
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
from matplotlib.patches import Polygon
import numpy as np


def beamline_get_ticks_locations_and_labels(o):
    ticks_locations = []
    ticks_labels = []
    for index, row in o.iterrows():
        if not re.match("DRIFT", index) \
            and not re.match(r".*\$", index) \
            and not re.match(r".*ENT",index) \
            and not re.match(".*EXT", index) \
            and not re.match(".*_.*", index):
                ticks_locations.append(row['S'])
                ticks_labels.append(str(index).split('_', 1)[0])
    return [ticks_locations, ticks_labels]


def draw_slits_and_markers(ax, context):
    pass


def twiss_plot(ax, twiss, context, with_plot=True, **kwargs):
    global palette
    ax1 = ax

    ## Plot Content
    if with_plot:
        filled_plot(ax1, twiss['S'], 0, 1000 * twiss['XMAXMONO'], palette['X'], True, alpha=0.8)
        filled_plot(ax1, twiss['S'], 0, 1000 * twiss['XMAX'], palette['X'], True, alpha=0.4)
        filled_plot(ax1, twiss['S'], 0, 2 * 1000 * twiss['XMAX'], palette['X'], True, alpha=0.2)

        filled_plot(ax1, twiss['S'], 0, 1000 * twiss['YMAX'], palette['Y'], True, alpha=0.4)
        filled_plot(ax1, twiss['S'], 0, 2 * 1000 * twiss['YMAX'], palette['Y'], True, alpha=0.2)

    ## Axes formatting
    ticks_locations, ticks_labels = beamline_get_ticks_locations_and_labels(twiss)
    ax1.tick_params(axis='both', which='major')
    ax1.tick_params(axis='x', labelsize=6)
    ax1.xaxis.set_major_locator(FixedLocator(ticks_locations))

    ax1.set_xlim([ticks_locations[0], ticks_locations[-1]])
    ax1.get_xaxis().set_tick_params(direction='out')
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=-45)
    ax1.yaxis.set_major_locator(MultipleLocator(10))
    ax1.yaxis.set_minor_locator(MultipleLocator(5))
    ax1.set_ylim([-45, 45])
    ax1.set_xlabel('s (m)')
    if kwargs.get("size_label", False):
        ax1.set_ylabel("{} beam size (mm)".format(kwargs['size_label']))
    else:
        ax1.set_ylabel(r'Beam size (mm)')
    ax1.grid(True, alpha=0.25)

    ax2 = ax1.twiny()
    ax2.set_xlim([ticks_locations[0], ticks_locations[-1]])
    ax2.get_xaxis().set_tick_params(direction='out')
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=-90)
    ax2.tick_params(axis='both', which='major')
    ax2.tick_params(axis='x', labelsize=6)
    ax2.xaxis.set_major_formatter(FixedFormatter(ticks_labels))
    ax2.xaxis.set_major_locator(FixedLocator(ticks_locations))

    ## Markers and slits
    draw_slits_and_markers(ax1, context)

    ## Arrows
    if kwargs.get("size_arrows", False):
        ax1.set_yticklabels([str(abs(x)) for x in ax1.get_yticks()])
        ax1.annotate('', xy=(-0.103, 0.97), xytext=(-0.103, 0.75),
                 arrowprops=dict(arrowstyle="->", color='k'), xycoords=ax1.transAxes)
        ax1.annotate('', xy=(-0.103, 0.25), xycoords='axes fraction', xytext=(-0.103, 0.03),
                 arrowprops=dict(arrowstyle="<-", color='k'))
        ax1.text(-0.126, 0.86, "Vertical", fontsize=7, rotation=90, transform=ax1.transAxes)
        ax1.text(-0.126, 0.22, "Horizontal", fontsize=7, rotation=90, transform=ax1.transAxes)


def draw_losses(ax, transmission, twiss):
    global palette
    ticks_locations, ticks_labels = beamline_get_ticks_locations_and_labels(twiss)
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


def filled_plot(ax, x, y0, y, c, fill=False, **kwargs):
    ax.plot(x, y, '.', markersize=0,
            markerfacecolor=c, markeredgecolor=c, color=c, **kwargs)
    if fill:
        ax.fill_between(x, y0, y, facecolor=c, linewidth=0.0, edgecolor=c, **kwargs)
