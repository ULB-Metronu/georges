import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib import patches
from numpy import linspace
from scipy.optimize import curve_fit


plt.ion()


# Define a function to make the ellipses
def ellipse(ra, rb, ang, x0, y0):

    theta = np.arange(0.0, 360.0, 1.0) * np.pi / 180.0

    xcenter, ycenter = x0, y0

    width = 2*ra
    height = 2*rb

    x = width * np.cos(theta)
    y = height * np.sin(theta)

    rtheta = np.radians(ang)
    rotation_matrix = np.array([
        [np.cos(rtheta), -np.sin(rtheta)],
        [np.sin(rtheta), np.cos(rtheta)],
    ])

    x, y = np.dot(rotation_matrix, np.array([x, y]))
    x += xcenter
    y += ycenter

    e1 = patches.Ellipse((xcenter, ycenter), width, height,
                         angle=ang,
                         linewidth=1,
                         fill=False,
                         linestyle='--',
                         edgecolor='black')

    return e1


def gauss(x, *p):
    a, mu, sigma = p
    return a*np.exp(-(x-mu)**2/(2.*sigma**2))


def draw2d_histo(fig, data, twiss_parameter, **kwargs):

    draw_ellipse = kwargs.get('draw_ellipse', False)

    # Define the x and y data
    x = data[data.columns[0]]
    y = data[data.columns[1]]

    xlims = [min(x), max(x)]
    ylims = [min(y), max(y)]

    # Set up default x and y limits
    if kwargs.get('xlim', [min(x), max(x)]):
        xlims = kwargs.get('xlim')

    if kwargs.get('ylim', [min(y), max(y)]):
        ylims = kwargs.get('ylim')
  
    # Define the locations for the axes
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.05
 
    # Set up the geometry of the three plots
    rect_beam = [left, bottom, width, height]  # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.25]  # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.25, height]  # dimensions of y-histogram
    rect_tab = [left_h, bottom_h, 0.25, 0.25]  # dimensions of tab

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
                          aspect='auto', cmap=kwargs.get('cmap', 'gray')))

    # Plot the beam plot contours
    contourcolor = 'black'
    xcenter = np.mean(x)
    ycenter = np.mean(y)
    ra = np.std(x)
    rb = np.std(y)
    ang = twiss_parameter['PHI']

    if draw_ellipse:

        ax_global.add_patch(ellipse(ra, rb, ang, xcenter, ycenter))
        ax_global.add_patch(ellipse(2*ra, 2*rb, ang, xcenter, ycenter))

        r1 = (xcenter+ra+0.25*ra)
        x1 = r1 * np.cos(np.deg2rad(ang))
        y1 = r1 * np.sin(np.deg2rad(ang))

        r2 = (xcenter+2*ra+0.25*ra)
        x2 = r2 * np.cos(np.deg2rad(ang))
        y2 = r2 * np.sin(np.deg2rad(ang))

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
    xbins = np.arange(xmin, xmax, (xmax-xmin)/nxbins)
    ybins = np.arange(ymin, ymax, (ymax-ymin)/nybins)

    # Plot the histograms
    nx, binsx, _ = ax_histx.hist(x, bins=xbins, color='blue', histtype='step', normed=True)
    ny, binsy, _ = ax_histy.hist(y, bins=ybins, orientation='horizontal', color='red', histtype='step', normed=True)

    # Set up the histogram limits
    ax_global.set_xlim(xlims)
    ax_global.set_ylim(ylims)
    ax_histx.set_xlim(xlims)
    ax_histy.set_ylim(ylims)

    if kwargs.get('gaussian_fit', False):

        bincenter = (binsx[:-1] + binsx[1:])/2
        coeff_x, var_matrix_x = curve_fit(gauss, bincenter, nx, p0=[1., x.mean(), x.std()])
        hist_fit = gauss(bincenter, *coeff_x)
        ax_histx.plot(bincenter, hist_fit, 'k--', linewidth=1)

        bincenter = (binsy[:-1] + binsy[1:]) / 2
        coeff_y, var_matrix_y = curve_fit(gauss, bincenter, ny, p0=[1., y.mean(), y.std()])
        hist_fit = gauss(bincenter, *coeff_y)
        ax_histy.plot(hist_fit, bincenter, 'k--', linewidth=1)

        error_x = np.sqrt(np.diag(var_matrix_x))
        error_y = np.sqrt(np.diag(var_matrix_y))

        print(f"GAUSSIAN MEAN X: {coeff_x[1]:{10}.{3}}+-{error_x[1]:{10}.{3}}\n"
              f"GAUSSIAN STD X: {coeff_x[2]:{10}.{3}}+- {error_x[2]:{10}.{3}}\n"
              f"GAUSSIAN MEAN Y: {coeff_y[1]:{10}.{3}}+- {error_y[1]:{10}.{3}} \n"
              f"GAUSSIAN STD Y: {coeff_y[2]:{10}.{3}}+- {error_y[2]:{10}.{3}}\n"
              )

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
        ax_tab.axhline(y0+0.4, color='k', linewidth=1)
        ax_tab.axhline(y0+0.8, color='k', linewidth=1)
        ax_tab.axhline(y0+0.83, color='k', linewidth=1)

        ax_tab.annotate('mean', xy=(x0+0.375, y0+0.9), xytext=(x0+0.375, y0+0.9),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate('std', xy=(x0+0.625, y0+0.9), xytext=(x0+0.625, y0+0.9),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate('median', xy=(x0+0.875, y0+0.9), xytext=(x0+0.875, y0+0.9),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)

        xname = '0'
        yname = '1'

        if kwargs.get('tab_names'):

            xname = kwargs.get('tab_names')[0]
            yname = kwargs.get('tab_names')[1]

        ax_tab.annotate(xname, xy=(x0+0.125, 0.6), xytext=(x0+0.125, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)

        ax_tab.annotate(yname, xy=(x0+0.125, 0.2), xytext=(x0+0.125, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)

        ax_tab.annotate(mean_x, xy=(x0+0.375, 0.6), xytext=(x0+0.375, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate(std_x, xy=(x0+0.625, 0.6), xytext=(x0+0.625, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate(median_x, xy=(x0+0.875, 0.6), xytext=(x0+0.875, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)

        ax_tab.annotate(mean_y, xy=(x0+0.375, 0.2), xytext=(x0+0.375, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate(std_y, xy=(x0+0.625, 0.2), xytext=(x0+0.625, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate(median_y, xy=(x0+0.875, 0.2), xytext=(x0+0.875, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
