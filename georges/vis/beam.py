import numpy as np
from numpy.linalg import eig
from numpy import linspace
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.ticker import NullFormatter
#from georges import statistics


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
    width = ra
    height = rb
    return patches.Ellipse((x0, y0), width, height,
                           angle=angle,
                           linewidth=kwargs.get('linewidth', 2),
                           fill=kwargs.get('fill', False),
                           linestyle=kwargs.get('linestyle', '--'),
                           edgecolor=kwargs.get('color', 'red'),
                           label=kwargs.get('label'),
                           )


def rotation_angle(eval, evec):
    if eval[0] * evec[0, 0] > eval[1] * evec[0, 1]:
        return np.degrees(np.arctan(evec[1, 0] / evec[0, 0]))
    else:
        return np.degrees(np.arctan(evec[1, 1] / evec[0, 1]))


def phase_space(fig, data, **kwargs):
    # SOURCE at
    # http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
    # Define the x and y data
    if data.columns[0] == 'X' or data.columns[0] == 'Y':
        unit_col_0 = '[mm]'
    else:
        unit_col_0 = '[mrad]'
    if data.columns[1] == 'X' or data.columns[1] == 'Y':
        unit_col_1 = '[mm]'
    else:
        unit_col_1 = '[mrad]'

    x = data[data.columns[0]]
    y = data[data.columns[1]]
    X_mean = np.mean(x)
    Y_mean = np.mean(y)

    xlims = [X_mean - 1.2 * np.max(np.abs(x)), X_mean + 1.2 * np.max(np.abs(x))]
    ylims = [Y_mean - 1.2 * np.max(np.abs(y)), Y_mean + 1.2 * np.max(np.abs(y))]

    CV = np.cov(x, y)
    [Eval, Evec] = eig(CV)
    chisq = [2.278868566, 5.991464547, 11.61828598]
    X_r, Y_r = [], []
    for i in range(0, len(chisq)):
        X_r.append(2 * (chisq[i] * Eval[0]) ** 0.5)
        Y_r.append(2 * (chisq[i] * Eval[1]) ** 0.5)

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
    nxbins = 100
    nybins = 100
    xbins = linspace(start=xmin, stop=xmax, num=nxbins)
    ybins = linspace(start=ymin, stop=ymax, num=nybins)
    h, xedges, yedges = np.histogram2d(y, x, bins=(ybins, xbins))
    # Plot the beam data
    _ = (ax_global.imshow(h, extent=[xmin, xmax, ymin, ymax],
                          interpolation='nearest', origin='lower',
                          aspect='auto', cmap=kwargs.get('cmap', 'gist_gray_r')))
    # Plot the beam plot contours
    contourcolor = 'black'

    ang = rotation_angle(Eval, Evec)
    ang_annotx = np.cos(np.radians(ang))
    ang_annoty = np.sin(np.radians(ang))
    if kwargs.get('draw_ellipse', True):
        ax_global.add_patch(ellipse(X_r[0], Y_r[0], ang, X_mean, Y_mean))
        ax_global.add_patch(ellipse(X_r[1], Y_r[1], ang, X_mean, Y_mean))
        ax_global.add_patch(ellipse(X_r[2], Y_r[2], ang, X_mean, Y_mean))
        ax_global.annotate('$1\\sigma$', xy=(X_r[0] / 2 * ang_annotx, Y_r[0] / 2 * ang_annoty),
                           xytext=(X_r[0] / 2 * ang_annotx, Y_r[0] / 2 * ang_annoty),
                           horizontalalignment='center', verticalalignment='center',
                           fontsize=15, color=contourcolor)
        ax_global.annotate('$2\\sigma$', xy=(X_r[1] / 2 * ang_annotx, Y_r[1] / 2 * ang_annoty),
                           xytext=(X_r[1] / 2 * ang_annotx, Y_r[1] / 2 * ang_annoty),
                           horizontalalignment='center', verticalalignment='center',
                           fontsize=15, color=contourcolor)
        ax_global.annotate('$3\\sigma$', xy=(X_r[2] / 2 * ang_annotx, Y_r[2] / 2 * ang_annoty),
                           xytext=(X_r[2] / 2 * ang_annotx, Y_r[2] / 2 * ang_annoty),
                           horizontalalignment='center', verticalalignment='center',
                           fontsize=15, color=contourcolor)
    # Plot the axes labels
    ax_global.set_xlabel(data.columns[0] + ' ' + unit_col_0, fontsize=18)
    ax_global.set_ylabel(data.columns[1] + ' ' + unit_col_1, fontsize=18)
    # Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax - xmin) / nxbins)
    ybins = np.arange(ymin, ymax, (ymax - ymin) / nybins)
    # Plot the histograms
    nx, binsx, _ = ax_histx.hist(x, bins=xbins, color='blue', histtype='step')
    ny, binsy, _ = ax_histy.hist(y, bins=ybins, orientation='horizontal', color='red', histtype='step')
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

        xname = data.columns[0] + '\n' + unit_col_0
        yname = data.columns[1] + '\n' + unit_col_1
        ax_tab.annotate(xname, xy=(x0 + 0.125, 0.6), xytext=(x0 + 0.125, 0.6),
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=11)
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


def phase_space_d(ax_global, ax_histx, ax_histy, ax_tab, data, elt, dim):
    # SOURCE at
    # http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
    # Define the x and y data

    if dim[0] == 'X' or dim[0] == 'Y':
        unit_col_0 = '[mm]'
    else:
        unit_col_0 = '[mrad]'
    if dim[1] == 'X' or dim[1] == 'Y':
        unit_col_1 = '[mm]'
    else:
        unit_col_1 = '[mrad]'

    #x = 1e3 * data[elt].distribution[dim[0]]
    #y = 1e3 * data[elt].distribution[dim[1]]

    x = 1e3 * data['BEAM_OUT'][elt][:,0]
    y = 1e3 * data['BEAM_OUT'][elt][:,2]


    X_mean = np.mean(x)
    Y_mean = np.mean(y)

    xlims = [X_mean - np.max(np.abs(x)), X_mean + np.max(np.abs(x))]
    ylims = [Y_mean - np.max(np.abs(y)), Y_mean + np.max(np.abs(y))]

    CV = np.cov(x, y)
    [Eval, Evec] = eig(CV)
    chisq = [2.278868566, 5.991464547, 11.61828598]
    X_r, Y_r = [], []
    for i in range(0, len(chisq)):
        X_r.append(2 * (chisq[i] * Eval[0]) ** 0.5)
        Y_r.append(2 * (chisq[i] * Eval[1]) ** 0.5)

    # Find the min/max of the data
    xmin, xmax, ymin, ymax = min(xlims), max(xlims), min(ylims), max(ylims)

    # Make the 'main' beam plot
    # Define the number of bins
    nxbins = 100
    nybins = 100
    xbins = linspace(start=xmin, stop=xmax, num=nxbins)
    ybins = linspace(start=ymin, stop=ymax, num=nybins)
    h, xedges, yedges = np.histogram2d(y, x, bins=(ybins, xbins))
    # Plot the beam data
    _ = (ax_global.imshow(h, extent=[xmin, xmax, ymin, ymax],
                          interpolation='nearest', origin='lower',
                          aspect='auto', cmap='gist_gray_r'))
    # Plot the beam plot contours
    ang = rotation_angle(Eval, Evec)
    ang_annotx = np.cos(np.radians(ang));
    ang_annoty = np.sin(np.radians(ang))
    ax_global.add_patch(ellipse(X_r[0], Y_r[0], ang, X_mean, Y_mean, color='red', label='$1\\sigma$'))
    ax_global.add_patch(ellipse(X_r[1], Y_r[1], ang, X_mean, Y_mean, color='blue', label='$2\\sigma$'))
    ax_global.add_patch(ellipse(X_r[2], Y_r[2], ang, X_mean, Y_mean, color='green', label='$3\\sigma$'))
    ax_global.legend()

    # Plot the axes labels
    ax_global.set_xlabel(dim[0] + ' ' + unit_col_0, fontsize=18)
    ax_global.set_ylabel(dim[1] + ' ' + unit_col_1, fontsize=18)
    # Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax - xmin) / nxbins)
    ybins = np.arange(ymin, ymax, (ymax - ymin) / nybins)
    # Plot the histograms
    nx, binsx, _ = ax_histx.hist(x, bins=xbins, color='blue', histtype='step')
    ny, binsy, _ = ax_histy.hist(y, bins=ybins, orientation='horizontal', color='red', histtype='step')
    # Set up the histogram limits
    ax_global.set_xlim(xlims)
    ax_global.set_ylim(ylims)
    ax_histx.set_xlim(xlims)
    ax_histy.set_ylim(ylims)

    limx = np.arange(xlims[0], xlims[1], (xlims[1] - xlims[0]) / nxbins)
    bin_centerx, fitresults_x = statistics.histogram_fit(x, bounds_binning=limx, verbose=False)
    ax_histx.plot(bin_centerx, fitresults_x.best_fit, 'k--', linewidth=1)
    limy = np.arange(ylims[0], ylims[1], (ylims[1] - ylims[0]) / nybins)
    bin_centery, fitresults_y = statistics.histogram_fit(y, bounds_binning=limy, verbose=False)
    ax_histy.plot(fitresults_y.best_fit, bin_centery, 'k--', linewidth=1)

    mean_x = str(round(x.mean(), 3))
    std_x = str(round(x.std(), 3))
    median_x = str(round(x.median(), 3))
    mean_y = str(round(y.mean(), 3))
    std_y = str(round(y.std(), 3))
    median_y = str(round(y.median(), 3))

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

    xname = dim[0] + '\n' + unit_col_0
    yname = dim[1] + '\n' + unit_col_1
    ax_tab.annotate(xname, xy=(x0 + 0.125, 0.6), xytext=(x0 + 0.125, 0.6),
                    horizontalalignment='center', verticalalignment='center',
                    fontsize=11)
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


def five_spot_map(bl_track0, bl_track1, bl_track2, bl_track3, bl_track4, figsize=(8, 8)):
    fig = plt.figure(figsize=figsize)
    try:
        # Order the 5 simulations outputs
        center = pd.DataFrame(bl_track0.line['BEAM']['ISO'].distribution['X'])
        blcorner = pd.DataFrame(bl_track1.line['BEAM']['ISO'].distribution['X'])
        brcorner = pd.DataFrame(bl_track2.line['BEAM']['ISO'].distribution['X'])
        tlcorner = pd.DataFrame(bl_track3.line['BEAM']['ISO'].distribution['X'])
        trcorner = pd.DataFrame(bl_track4.line['BEAM']['ISO'].distribution['X'])
        dataX = pd.concat([center, blcorner, brcorner, tlcorner, trcorner], ignore_index=True)

        center = pd.DataFrame(bl_track0.line['BEAM']['ISO'].distribution['Y'])
        blcorner = pd.DataFrame(bl_track1.line['BEAM']['ISO'].distribution['Y'])
        brcorner = pd.DataFrame(bl_track2.line['BEAM']['ISO'].distribution['Y'])
        tlcorner = pd.DataFrame(bl_track3.line['BEAM']['ISO'].distribution['Y'])
        trcorner = pd.DataFrame(bl_track4.line['BEAM']['ISO'].distribution['Y'])
        dataY = pd.concat([center, blcorner, brcorner, tlcorner, trcorner], ignore_index=True)
    except:
        print('Error, one simulation is missing. You have to put the 5 simulations outputs in the functions')

    # Create the 2D histogram with the 5 spots.
    _ = plt.hist2d(dataX['X'].values, dataY['Y'].values, bins=400, cmap='gist_gray_r')
    data_X = []
    data_X.append(1e3 * bl_track0.line['BEAM']['ISO'].std['X'])
    data_X.append(1e3 * bl_track1.line['BEAM']['ISO'].std['X'])
    data_X.append(1e3 * bl_track2.line['BEAM']['ISO'].std['X'])
    data_X.append(1e3 * bl_track3.line['BEAM']['ISO'].std['X'])
    data_X.append(1e3 * bl_track4.line['BEAM']['ISO'].std['X'])

    data_S = []
    data_S.append(100 * np.abs(bl_track0.line['BEAM']['ISO'].std['X'] - bl_track0.line['BEAM']['ISO'].std['Y']) / (
    bl_track0.line['BEAM']['ISO'].std['X'] + bl_track0.line['BEAM']['ISO'].std['Y']))
    data_S.append(100 * np.abs(bl_track1.line['BEAM']['ISO'].std['X'] - bl_track1.line['BEAM']['ISO'].std['Y']) / (
    bl_track1.line['BEAM']['ISO'].std['X'] + bl_track1.line['BEAM']['ISO'].std['Y']))
    data_S.append(100 * np.abs(bl_track2.line['BEAM']['ISO'].std['X'] - bl_track2.line['BEAM']['ISO'].std['Y']) / (
    bl_track2.line['BEAM']['ISO'].std['X'] + bl_track2.line['BEAM']['ISO'].std['Y']))
    data_S.append(100 * np.abs(bl_track3.line['BEAM']['ISO'].std['X'] - bl_track3.line['BEAM']['ISO'].std['Y']) / (
    bl_track3.line['BEAM']['ISO'].std['X'] + bl_track3.line['BEAM']['ISO'].std['Y']))
    data_S.append(100 * np.abs(bl_track4.line['BEAM']['ISO'].std['X'] - bl_track4.line['BEAM']['ISO'].std['Y']) / (
    bl_track4.line['BEAM']['ISO'].std['X'] + bl_track4.line['BEAM']['ISO'].std['Y']))

    data_T = []
    data_T.append(np.degrees(ellipse_angle_of_rotation(fitEllipse(1e3 * bl_track0.line['BEAM']['ISO'].distribution['X'],
                                                                  1e3 * bl_track0.line['BEAM']['ISO'].distribution[
                                                                      'Y']))))
    data_T.append(np.degrees(ellipse_angle_of_rotation(fitEllipse(1e3 * bl_track1.line['BEAM']['ISO'].distribution['X'],
                                                                  1e3 * bl_track1.line['BEAM']['ISO'].distribution[
                                                                      'Y']))))
    data_T.append(np.degrees(ellipse_angle_of_rotation(fitEllipse(1e3 * bl_track2.line['BEAM']['ISO'].distribution['X'],
                                                                  1e3 * bl_track2.line['BEAM']['ISO'].distribution[
                                                                      'Y']))))
    data_T.append(np.degrees(ellipse_angle_of_rotation(fitEllipse(1e3 * bl_track3.line['BEAM']['ISO'].distribution['X'],
                                                                  1e3 * bl_track3.line['BEAM']['ISO'].distribution[
                                                                      'Y']))))
    data_T.append(np.degrees(ellipse_angle_of_rotation(fitEllipse(1e3 * bl_track4.line['BEAM']['ISO'].distribution['X'],
                                                                  1e3 * bl_track4.line['BEAM']['ISO'].distribution[
                                                                      'Y']))))

    plt.annotate('$\sigma_X = {}$mm \n $S = {}$% \n $\phi={} °$'.format(round(data_X[1], 2), round(data_S[1], 2),
                                                                        round(data_T[1], 2))
                 , xy=(- 0.1 + 0.03, 0.1), xytext=(- 0.1 + 0.03, 0.1),
                 horizontalalignment='center', verticalalignment='center',
                 fontsize=12)
    plt.annotate('$\sigma_X = {}$mm \n $S = {}$% \n $\phi={} °$'.format(round(data_X[4], 2), round(data_S[4], 2),
                                                                        round(data_T[4], 2))
                 , xy=(0.1 - 0.03, 0.1), xytext=(0.1 - 0.03, 0.1),
                 horizontalalignment='center', verticalalignment='center',
                 fontsize=12)
    plt.annotate('$\sigma_X = {}$mm \n $S = {}$% \n $\phi={} °$'.format(round(data_X[2], 2), round(data_S[2], 2),
                                                                        round(data_T[2], 2))
                 , xy=(- 0.1 + 0.03, -0.1), xytext=(- 0.1 + 0.03, -0.1),
                 horizontalalignment='center', verticalalignment='center',
                 fontsize=12)
    plt.annotate('$\sigma_X = {}$mm \n $S = {}$% \n $\phi={} °$'.format(round(data_X[3], 2), round(data_S[3], 2),
                                                                        round(data_T[3], 2))
                 , xy=(0.1 - 0.03, -0.1), xytext=(0.1 - 0.03, -0.1),
                 horizontalalignment='center', verticalalignment='center',
                 fontsize=12)
    plt.annotate('$\sigma_X = {}$mm \n $S = {}$% \n $\phi={} °$'.format(round(data_X[0], 2), round(data_S[0], 2),
                                                                        round(data_T[0], 2))
                 , xy=(0.0 + 0.03, 0.0), xytext=(0.0 + 0.03, 0.0),
                 horizontalalignment='center', verticalalignment='center',
                 fontsize=12)