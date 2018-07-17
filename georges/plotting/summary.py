import matplotlib.pyplot as plt
from matplotlib.ticker import *
from .beam import phase_space_d
from .common import prepare
from .aperture import aperture
from .tracking import tracking
from .losses import losses


def summary(bl, bl_track, context, element='ISO'):
    fig = plt.figure(figsize=(20, 20))
    # Define the left bottom corner block
    space = 0.02
    width, height = 0.25, 0.25
    left, hwidth = 0.05, width / 2
    bottom, hheight = 0.05, height / 2
    bottom_h = left_h = left + width + space
    # Set up the geometry of the three plots
    rect_beam = [left, bottom, width, height]  # dimensions of temp plot
    rect_histx = [left, bottom_h, width, hheight]  # dimensions of x-histogram
    rect_histy = [left_h, bottom, hwidth, height]  # dimensions of y-histogram
    rect_tab = [left_h, bottom_h, hwidth, hheight]  # dimensions of tab
    # Make the three plots
    ax_global = fig.add_axes(rect_beam)  # beam plot
    ax_histx = fig.add_axes(rect_histx)  # x histogram
    ax_histy = fig.add_axes(rect_histy)  # y histogram
    ax_histx.set_ylabel("Counts")
    ax_histx.grid(True)
    ax_histy.set_xlabel("Counts")
    ax_histy.grid(True)
    nullfmt = NullFormatter()
    ax_histx.xaxis.set_major_formatter(nullfmt)
    ax_histy.yaxis.set_major_formatter(nullfmt)
    ax_tab = fig.add_axes(rect_tab)  # y histogram
    ax_tab.tick_params(labelbottom='off', labelleft='off', left='off', bottom='off')

    dim = ['Y', 'PY']
    phase_space_d(ax_global, ax_histx, ax_histy, ax_tab, bl_track.line['BEAM'], element, dim)

    # Define the right bottom corner block
    x_bound = left + width + space + hwidth + 3 * space
    n_left = x_bound
    left_h = n_left + width + space
    bottom_h = bottom + height + space
    # Set up the geometry of the three plots
    rect_beam = [n_left, bottom, width, height]  # dimensions of temp plot
    rect_histx = [n_left, bottom_h, width, hheight]  # dimensions of x-histogram
    rect_histy = [left_h, bottom, hwidth, height]  # dimensions of y-histogram
    rect_tab = [left_h, bottom_h, hwidth, hheight]  # dimensions of tab
    # Make the three plots
    ax_global = fig.add_axes(rect_beam)  # beam plot
    ax_histx = fig.add_axes(rect_histx)  # x histogram
    ax_histy = fig.add_axes(rect_histy)  # y histogram
    ax_histx.set_ylabel("Counts")
    ax_histx.grid(True)
    ax_histy.set_xlabel("Counts")
    ax_histy.grid(True)
    nullfmt = NullFormatter()
    ax_histx.xaxis.set_major_formatter(nullfmt)
    ax_histy.yaxis.set_major_formatter(nullfmt)
    ax_tab = fig.add_axes(rect_tab)  # y histogram
    ax_tab.tick_params(labelbottom='off', labelleft='off', left='off', bottom='off')

    dim = ['X', 'PX']
    phase_space_d(ax_global, ax_histx, ax_histy, ax_tab, bl_track.line['BEAM'], element, dim)

    # Define the right top corner block
    y_bound = bottom + height + space + hheight + 3 * space
    n_bottom = y_bound
    left_h = n_left + width + space
    bottom_h = n_bottom + height + space
    # Set up the geometry of the three plots
    rect_beam = [n_left, n_bottom, width, height]  # dimensions of temp plot
    rect_histx = [n_left, bottom_h, width, hheight]  # dimensions of x-histogram
    rect_histy = [left_h, n_bottom, hwidth, height]  # dimensions of y-histogram
    rect_tab = [left_h, bottom_h, hwidth, hheight]  # dimensions of tab
    # Make the three plots
    ax_global = fig.add_axes(rect_beam)  # beam plot
    ax_histx = fig.add_axes(rect_histx)  # x histogram
    ax_histy = fig.add_axes(rect_histy)  # y histogram
    ax_histx.set_ylabel("Counts")
    ax_histx.grid(True)
    ax_histy.set_xlabel("Counts")
    ax_histy.grid(True)
    nullfmt = NullFormatter()
    ax_histx.xaxis.set_major_formatter(nullfmt)
    ax_histy.yaxis.set_major_formatter(nullfmt)
    ax_tab = fig.add_axes(rect_tab)  # y histogram

    dim = ['X', 'Y']
    phase_space_d(ax_global, ax_histx, ax_histy, ax_tab, bl_track.line['BEAM'], element, dim)

    # Define the left top corner block
    n_height = (height + space + hheight) / 3
    n_width = width + space + hwidth
    bottom_h = n_bottom - 1.5 * space + n_height + space
    bottom_hh = bottom_h + n_height + space
    # Set up the geometry of the three plots
    rect_beam_X = [left, bottom_hh, n_width, n_height]
    rect_beam_Y = [left, bottom_h, n_width, n_height]
    rect_trans = [left, n_bottom - 1.5 * space, n_width, n_height]
    # Make the three plots
    ax_beam_X = fig.add_axes(rect_beam_X)  # beam plot
    ax_beam_Y = fig.add_axes(rect_beam_Y)  # x histogram
    ax_trans = fig.add_axes(rect_trans)  # y histogram

    prepare(ax_beam_X, bl)
    aperture(ax_beam_X, bl, context=context, plane='X')
    tracking(ax_beam_X, bl_track, context=context, plane='X', halo=True, halo99=True, std=True, mean=True)

    prepare(ax_beam_Y, bl, print_label=False)
    aperture(ax_beam_Y, bl, context=context, plane='Y')
    tracking(ax_beam_Y, bl_track, context=context, plane='Y', halo=True, halo99=True, std=True, mean=True)

    prepare(ax_trans, bl, print_label=False)
    losses_ = losses(ax_trans, bl_track, log=False)

    ax_beam_X.set_ylabel("Horizontal beam size [mm]")
    ax_beam_Y.set_ylabel("Vertical beam size [mm]")
    ax_beam_X.set_xticklabels([])
    ax_beam_Y.set_xticklabels([])
    ax_beam_X.set_xlabel('')
    ax_beam_Y.set_xlabel('')
    ax_beam_X.set_ylim([-40, 40])
    ax_beam_Y.set_ylim([-40, 40])