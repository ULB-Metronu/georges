import numpy as np
from numpy import linspace
from matplotlib import patches
from matplotlib.ticker import NullFormatter
from .. import statistics


class BeamPlottingException(Exception):
    """Exception raised for errors in the beam plotting module."""

    def __init__(self, m):
        self.message = m


def ellipse(ra, rb, angle, x0, y0, **kwargs):
    """
    Create an ellipse from beam parameters.
    :param ra: semi-major axis
    :param rb: semi-minor axis
    :param angle: oritentation angle
    :param x0: center X coordinate
    :param y0: center Y coordinate
    :return: `matplotlib.patches.Ellipse` object
    """

    width = 2 * ra
    height = 2 * rb

    return patches.Ellipse((x0, y0), width, height,
                           angle=np.rad2deg(angle),
                           linewidth=kwargs.get('linewidth', 1),
                           fill=kwargs.get('fill', False),
                           linestyle=kwargs.get('linestyle', '--'),
                           edgecolor=kwargs.get('color', 'black'))


def phase_space(fig, data, **kwargs):

    # Define the x and y data
    x = data[data.columns[0]]
    y = data[data.columns[1]]

    xlims = [min(x), max(x)]
    ylims = [min(y), max(y)]

    if kwargs.get('twiss_angle') is None:
        raise BeamPlottingException("Error no angle specified.")

    # Set up default x and y limits
    if kwargs.get('xlims') is not None:
        xlims = kwargs.get('xlims')

    if kwargs.get('ylims') is not None:
        ylims = kwargs.get('ylims')

    # Define the locations for the axes
    left, width = 0.12, 0.6
    bottom, height = 0.12, 0.6
    bottom_h = left_h = left + width + 0.05

    # Set up the geometry of the three plots
    rect_beam = [left, bottom, width, height]  # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.25]  # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.4, height]  # dimensions of y-histogram
    rect_tab = [left_h, bottom_h, 0.4, 0.25]  # dimensions of tab

    # Make the three plots
    ax_global = fig.add_axes(rect_beam)  # beam plot
    ax_histx = fig.add_axes(rect_histx)  # x histogram
    ax_histy = fig.add_axes(rect_histy)  # y histogram

    ax_histx.set_ylabel("Counts")
    ax_histx.grid(True)

    ax_histy.set_xlabel("Counts")
    ax_histy.grid(True)

    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    ax_histx.xaxis.set_major_formatter(nullfmt)
    ax_histy.yaxis.set_major_formatter(nullfmt)

    # Find the min/max of the data
    xmin = min(xlims)
    xmax = max(xlims)
    ymin = min(ylims)
    ymax = max(ylims)

    # Make the 'main' beam plot
    # Define the number of bins
    if kwargs.get('nbins'):
        nxbins = kwargs.get('nbins')[0]
        nybins = kwargs.get('nbins')[1]
    else:
        nxbins = 50
        nybins = 50

    xbins = linspace(start=xmin, stop=xmax, num=nxbins)
    ybins = linspace(start=ymin, stop=ymax, num=nybins)
    h, xedges, yedges = np.histogram2d(y, x, bins=(ybins, xbins))

    # Plot the beam data
    _ = (ax_global.imshow(h, extent=[xmin, xmax, ymin, ymax],
                          interpolation='nearest', origin='lower',
                          aspect='auto', cmap=kwargs.get('cmap', 'gist_gray_r')))

    # Plot the beam plot contours
    contourcolor = 'black'
    xcenter = np.mean(x)
    ycenter = np.mean(y)
    ra = np.std(x)
    rb = np.std(y)
    ang = kwargs.get('twiss_angle')

    if kwargs.get('draw_ellipse', False):
        ax_global.add_patch(ellipse(ra, rb, ang, xcenter, ycenter))
        ax_global.add_patch(ellipse(2 * ra, 2 * rb, ang, xcenter, ycenter))

        x1, y1 = kwargs.get("xtext_pos")
        x2, y2 = kwargs.get("ytext_pos")

        ax_global.annotate('$1\\sigma$', xy=(x1, y1), xytext=(x1, y1),
                           horizontalalignment='center', verticalalignment='center',
                           fontsize=15, color=contourcolor)

        ax_global.annotate('$2\\sigma$', xy=(x2, y2), xytext=(x2, y2),
                           horizontalalignment='center', verticalalignment='center',
                           fontsize=15, color=contourcolor)

    # Plot the axes labels

    ax_global.set_xlabel(data.columns[0], fontsize=25)
    ax_global.set_ylabel(data.columns[1], fontsize=25)

    # Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax - xmin) / nxbins)
    ybins = np.arange(ymin, ymax, (ymax - ymin) / nybins)

    # Plot the histograms
    nx, binsx, _ = ax_histx.hist(x, bins=xbins, color='blue', histtype='step', normed=False)
    ny, binsy, _ = ax_histy.hist(y, bins=ybins, orientation='horizontal', color='red', histtype='step', normed=False)

    # Set up the histogram limits
    ax_global.set_xlim(xlims)
    ax_global.set_ylim(ylims)
    ax_histx.set_xlim(xlims)
    ax_histy.set_ylim(ylims)

    if kwargs.get('gaussian_fit', False):

        limx = np.arange(xlims[0], xlims[1], (xlims[1] - xlims[0]) / nxbins)
        bin_centerx, fitresults_x = statistics.histogram_fit(x, bounds_binning=limx,
                                                             verbose=kwargs.get("verbose", False))
        ax_histx.plot(bin_centerx, fitresults_x.best_fit, 'k--', linewidth=1)

        limy = np.arange(ylims[0], ylims[1], (ylims[1] - ylims[0]) / nybins)
        bin_centery, fitresults_y = statistics.histogram_fit(y, bounds_binning=limy,
                                                             verbose=kwargs.get("verbose", False))
        ax_histy.plot(fitresults_y.best_fit, bin_centery, 'k--', linewidth=1)

    if kwargs.get('draw_tab', False):

        mean_x = str(round(x.mean(), 3))
        std_x = str(round(x.std(), 3))
        median_x = str(round(x.median(), 3))
        mean_y = str(round(y.mean(), 3))
        std_y = str(round(y.std(), 3))
        median_y = str(round(y.median(), 3))

        ax_tab = fig.add_axes(rect_tab)  # y histogram

        ax_tab.tick_params(labelbottom='off', labelleft='off', left='off', bottom='off')

        x0 = ax_tab.get_xlim()[0]
        y0 = ax_tab.get_ylim()[0]

        ax_tab.axvline(x0 + 0.25, color='k', linewidth=1)
        ax_tab.axvline(x0 + 0.5, color='k', linewidth=1)
        ax_tab.axvline(x0 + 0.75, color='k', linewidth=1)
        ax_tab.axhline(y0 + 0.4, color='k', linewidth=1)
        ax_tab.axhline(y0 + 0.8, color='k', linewidth=1)
        ax_tab.axhline(y0 + 0.83, color='k', linewidth=1)

        ax_tab.annotate('mean', xy=(x0 + 0.375, y0 + 0.9), xytext=(x0 + 0.375, y0 + 0.9),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=12)
        ax_tab.annotate('std ', xy=(x0 + 0.625, y0 + 0.9), xytext=(x0 + 0.625, y0 + 0.9),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=12)
        ax_tab.annotate('median', xy=(x0 + 0.875, y0 + 0.9), xytext=(x0 + 0.875, y0 + 0.9),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=12)

        xname = '0'
        yname = '1'

        if kwargs.get('tab_names'):
            xname = kwargs.get('tab_names')[0]
            yname = kwargs.get('tab_names')[1]

        ax_tab.annotate(xname, xy=(x0 + 0.125, 0.6), xytext=(x0 + 0.125, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=12)

        ax_tab.annotate(yname, xy=(x0 + 0.125, 0.2), xytext=(x0 + 0.125, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=12)

        ax_tab.annotate(mean_x, xy=(x0 + 0.375, 0.6), xytext=(x0 + 0.375, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=12)
        ax_tab.annotate(std_x, xy=(x0 + 0.625, 0.6), xytext=(x0 + 0.625, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=12)
        ax_tab.annotate(median_x, xy=(x0 + 0.875, 0.6), xytext=(x0 + 0.875, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=12)

        ax_tab.annotate(mean_y, xy=(x0 + 0.375, 0.2), xytext=(x0 + 0.375, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=12)
        ax_tab.annotate(std_y, xy=(x0 + 0.625, 0.2), xytext=(x0 + 0.625, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=12)
        ax_tab.annotate(median_y, xy=(x0 + 0.875, 0.2), xytext=(x0 + 0.875, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=12)
