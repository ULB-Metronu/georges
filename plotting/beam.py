import numpy as np
from numpy import linspace
from matplotlib import patches
from matplotlib.ticker import NullFormatter, MaxNLocator
from .. import statistics as stat


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
    theta = np.radians(np.arange(0.0, 360.0, 1.0))

    width = 2 * ra
    height = 2 * rb

    x = width * np.cos(theta)
    y = height * np.sin(theta)

    rtheta = np.radians(angle)
    rotation_matrix = np.array([
        [np.cos(rtheta), -np.sin(rtheta)],
        [np.sin(rtheta), np.cos(rtheta)],
    ])

    x, y = np.dot(rotation_matrix, np.array([x, y]))
    x += x0
    y += y0

    return patches.Ellipse((x0, y0), width, height,
                           angle=angle,
                           linewidth=kwargs.get('linewidth', 1),
                           fill=kwargs.get('fill', False),
                           linestyle=kwargs.get('linestyle', '--'),
                           edgecolor=kwargs.get('color', 'black'))


def draw2d_histo(fig, data, twiss_angle, **kwargs):

    # Define the x and y data
    x = data[data.columns[0]]
    y = data[data.columns[1]]

    xlims = [min(x), max(x)]
    ylims = [min(y), max(y)]

    # Set up default x and y limits
    if kwargs.get('xlim') is not None:
        xlims = kwargs.get('xlim')

    if kwargs.get('ylim') is not None:
        ylims = kwargs.get('ylim')

    # Define the locations for the axes
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left + width + 0.05

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
                          aspect='auto', cmap=kwargs.get('cmap', 'gist_gray_r')))

    # Plot the beam plot contours
    contourcolor = 'black'
    xcenter = np.mean(x)
    ycenter = np.mean(y)
    ra = np.std(x)
    rb = np.std(y)
    ang = np.rad2deg(twiss_angle)

    if kwargs.get('draw_ellipse', False):
        ax_global.add_patch(ellipse(ra, rb, ang, xcenter, ycenter))
        ax_global.add_patch(ellipse(2 * ra, 2 * rb, ang, xcenter, ycenter))

        r1 = (xcenter + ra + 0.3 * ra)
        x1 = r1 * np.cos(np.deg2rad(ang))
        y1 = r1 * np.sin(np.deg2rad(ang))

        r2 = (xcenter + 2 * ra + 0.3 * ra)
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
    xbins = np.arange(xmin, xmax, (xmax - xmin) / nxbins)
    ybins = np.arange(ymin, ymax, (ymax - ymin) / nybins)

    # Plot the histograms
    nx, binsx, _ = ax_histx.hist(x, bins=xbins, color='blue', histtype='step', normed=True)
    ny, binsy, _ = ax_histy.hist(y, bins=ybins, orientation='horizontal', color='red', histtype='step', normed=True)

    # Set up the histogram limits
    ax_global.set_xlim(xlims)
    ax_global.set_ylim(ylims)
    ax_histx.set_xlim(xlims)
    ax_histy.set_ylim(ylims)

    if kwargs.get('gaussian_fit', False):
        limx = np.arange(xlims[0], xlims[1], (xlims[1] - xlims[0]) / nxbins)
        bin_centerx, fitresults_x = stat.gaussian_fit(x, default_lim=False, lim=limx)
        ax_histx.plot(bin_centerx, fitresults_x.best_fit, 'k--', linewidth=1)

        limy = np.arange(ylims[0], ylims[1], (ylims[1] - ylims[0]) / nybins)
        bin_centery, fitresults_y = stat.gaussian_fit(y, default_lim=False, lim=limy)
        ax_histy.plot(fitresults_y.best_fit, bin_centery, 'k--', linewidth=1)

        gaussian_meanx = float(fitresults_x.params['cen'].value)
        gaussian_meany = float(fitresults_y.params['cen'].value)
        gaussian_stdx = float(fitresults_x.params['wid'].value)
        gaussian_stdy = float(fitresults_y.params['wid'].value)

        gaussian_errormeanx = float(fitresults_x.params['cen'].stderr)
        gaussian_errormeany = float(fitresults_y.params['cen'].stderr)
        gaussian_errorstdx = float(fitresults_x.params['wid'].stderr)
        gaussian_errorstdy = float(fitresults_y.params['wid'].stderr)

        relative_error_meanx = 100 * gaussian_errormeanx / gaussian_meanx
        relative_error_meany = 100 * gaussian_errormeany / gaussian_meany
        relative_error_stdx = 100 * gaussian_errorstdx / gaussian_stdx
        relative_error_stdy = 100 * gaussian_errorstdy / gaussian_stdy

        print(f"GAUSSIAN MEAN X: {gaussian_meanx:{10}.{3}}"
              f"+/- {gaussian_errormeanx:{10}.{3}}"
              f"({relative_error_meanx:{10}.{3}}%)\n"
              f"GAUSSIAN STD X: {gaussian_stdx:{10}.{3}}"
              f"+/- {gaussian_errorstdx:{10}.{3}}"
              f"({relative_error_stdx:{10}.{3}}%)\n"
              f"GAUSSIAN MEAN Y: {gaussian_meany:{10}.{3}}"
              f"+/- {gaussian_errormeany:{10}.{3}}"
              f"({relative_error_meany:{10}.{3}}%)\n"
              f"GAUSSIAN STD Y: {gaussian_stdy:{10}.{3}}"
              f"+/- {gaussian_errorstdy:{10}.{3}}"
              f"({relative_error_stdy:{10}.{3}}%)\n"
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
        ax_tab.axhline(y0 + 0.4, color='k', linewidth=1)
        ax_tab.axhline(y0 + 0.8, color='k', linewidth=1)
        ax_tab.axhline(y0 + 0.83, color='k', linewidth=1)

        ax_tab.annotate('mean', xy=(x0 + 0.375, y0 + 0.9), xytext=(x0 + 0.375, y0 + 0.9),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate('std', xy=(x0 + 0.625, y0 + 0.9), xytext=(x0 + 0.625, y0 + 0.9),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate('median', xy=(x0 + 0.875, y0 + 0.9), xytext=(x0 + 0.875, y0 + 0.9),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)

        xname = '0'
        yname = '1'

        if kwargs.get('tab_names'):
            xname = kwargs.get('tab_names')[0]
            yname = kwargs.get('tab_names')[1]

        ax_tab.annotate(xname, xy=(x0 + 0.125, 0.6), xytext=(x0 + 0.125, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)

        ax_tab.annotate(yname, xy=(x0 + 0.125, 0.2), xytext=(x0 + 0.125, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)

        ax_tab.annotate(mean_x, xy=(x0 + 0.375, 0.6), xytext=(x0 + 0.375, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate(std_x, xy=(x0 + 0.625, 0.6), xytext=(x0 + 0.625, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate(median_x, xy=(x0 + 0.875, 0.6), xytext=(x0 + 0.875, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)

        ax_tab.annotate(mean_y, xy=(x0 + 0.375, 0.2), xytext=(x0 + 0.375, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate(std_y, xy=(x0 + 0.625, 0.2), xytext=(x0 + 0.625, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)
        ax_tab.annotate(median_y, xy=(x0 + 0.875, 0.2), xytext=(x0 + 0.875, 0.2),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=6)


# Define a function to make the ellipses

# def ellipse(ra, rb, ang, x0, y0, Nb=100):
#     xpos, ypos = x0, y0
#     radm, radn = ra, rb
#     an = ang
#     co, si = np.cos(an), np.sin(an)
#     the = linspace(0, 2 * np.pi, Nb)
#     X = radm * np.cos(the) * co - si * radn * np.sin(the) + xpos
#     Y = radm * np.cos(the) * si + co * radn * np.sin(the) + ypos
#     return X, Y


def make2Dplot(fig, Data_BEAMX, Data_BEAMY, Nbinx, Nbiny):
    # Define the x and y data
    # For example just using random numbers

    x = Data_BEAMX
    y = Data_BEAMY

    # Set up default x and y limits
    xlims = [min(x), max(x)]
    ylims = [min(y), max(y)]

    # Define the locations for the axes
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left + width + 0.02

    # Set up the geometry of the three plots
    rect_beam = [left, bottom, width, height]  # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.25]  # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.25, height]  # dimensions of y-histogram

    # Make the three plots
    axBeam = fig.add_axes(rect_beam)  # beam plot
    mean_DataX = str(round(Data_BEAMX.mean(), 3))
    std_DataX = str(round(Data_BEAMX.std(), 3))
    mean_DataY = str(round(Data_BEAMY.mean(), 3))
    std_DataY = str(round(Data_BEAMY.std(), 3))

    axHistx = fig.add_axes(rect_histx)  # x histogram
    axHistx.set_ylabel("Counts")
    axHistx.set_title('Mean : ' + mean_DataX + ' std : ' + std_DataX)
    axHistx.grid(True)

    axHisty = fig.add_axes(rect_histy)  # y histogram
    axHisty.set_xlabel("Counts")
    axHisty.set_title('Mean : ' + mean_DataY + ' std : ' + std_DataY, rotation=270, x=1.08, y=0.75)
    axHisty.grid(True)

    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # Find the min/max of the data
    xmin = min(xlims)
    xmax = max(xlims)
    ymin = min(ylims)
    ymax = max(ylims)

    # Make the 'main' beam plot
    # Define the number of bins
    nxbins = Nbinx
    nybins = Nbiny

    xbins = linspace(start=xmin, stop=xmax, num=nxbins)
    ybins = linspace(start=ymin, stop=ymax, num=nybins)
    xcenter = (xbins[0:-1] + xbins[1:]) / 2.0
    ycenter = (ybins[0:-1] + ybins[1:]) / 2.0
    aspectratio = 1.0 * (xmax - 0) / (1.0 * ymax - 0)

    H, xedges, yedges = np.histogram2d(y, x, bins=(ybins, xbins))
    X = xcenter
    Y = ycenter
    Z = H

    # Plot the beam data
    cax = (axBeam.imshow(H, extent=[xmin, xmax, ymin, ymax],
                         interpolation='nearest', origin='lower', aspect='auto'))

    print('xmin : ' + str(xmin) + ' xmax : ' + str(xmax) + ' ymin: ' + str(ymin) + ' ymax: ' + str(ymax))
    # Plot the beam plot contours
    contourcolor = 'white'
    xcenter = np.mean(x)
    ycenter = np.mean(y)
    ra = np.std(x)
    rb = np.std(y)
    ang = 0  ##To change for rotated ellipse : call georges.phys

    X, Y = ellipse(ra, rb, ang, xcenter, ycenter)
    axBeam.plot(X, Y, "k:", ms=1, linewidth=2.0)
    axBeam.annotate('$1\\sigma$', xy=(X[15], Y[15]), xycoords='data', xytext=(10, 10),
                    textcoords='offset points', horizontalalignment='right',
                    verticalalignment='bottom', fontsize=25)

    # X,Y=ellipse(2*ra,2*rb,ang,xcenter,ycenter)
    # axBeam.plot(X,Y,"k:",color = contourcolor,ms=1,linewidth=2.0)
    # axBeam.annotate('$2\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
    #             textcoords='offset points',horizontalalignment='right',
    #             verticalalignment='bottom',fontsize=25, color = contourcolor)

    X, Y = ellipse(3 * ra, 3 * rb, ang, xcenter, ycenter)
    axBeam.plot(X, Y, "k:", color=contourcolor, ms=1, linewidth=2.0)
    axBeam.annotate('$3\\sigma$', xy=(X[15], Y[15]), xycoords='data', xytext=(10, 10),
                    textcoords='offset points', horizontalalignment='right',
                    verticalalignment='bottom', fontsize=25, color=contourcolor)

    # Plot the axes labels
    axBeam.set_xlabel(Data_BEAMX.name, fontsize=25)
    axBeam.set_ylabel(Data_BEAMY.name, fontsize=25)

    # Make the tickmarks pretty
    ticklabels = axBeam.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(18)
        label.set_family('serif')

    ticklabels = axBeam.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(18)
        label.set_family('serif')

    # Set up the plot limits
    axBeam.set_xlim(xlims)
    axBeam.set_ylim(ylims)

    # Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax - xmin) / nxbins)
    ybins = np.arange(ymin, ymax, (ymax - ymin) / nybins)

    # Plot the histograms
    axHistx.hist(x, bins=xbins, color='blue', histtype='step', normed=True)
    axHisty.hist(y, bins=ybins, orientation='horizontal', color='red', histtype='step', normed=True)
    # Set up the histogram limits
    axHistx.set_xlim(min(x), max(x))
    axHisty.set_ylim(min(y), max(y))

    # Make the tickmarks pretty
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family('serif')

    # Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family('serif')

    # Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))
