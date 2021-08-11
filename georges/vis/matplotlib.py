"""Matplotlib plotting module for Manzoni.

TODO
"""
from __future__ import annotations
import numpy as _np
from numpy.linalg import eig
from numpy import linspace
import pandas as _pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import *
from matplotlib.ticker import NullFormatter

from lmfit.models import GaussianModel

from georges_core.vis import MatplotlibArtist as _MatplotlibArtist

PALETTE = {
    'solarized': {
        'base03': '#002b36',
        'base02': '#073642',
        'base01': '#586e75',
        'base00': '#657b83',
        'base0': '#839496',
        'base1': '#93a1a1',
        'base2': '#eee8d5',
        'base3': '#fdf6e3',
        'yellow': '#b58900',
        'orange': '#cb4b16',
        'red': '#dc322f',
        'magenta': '#d33682',
        'violet': '#6c71c4',
        'blue': '#268bd2',
        'cyan': '#2aa198',
        'green': '#859900'
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
palette['Quadrupole'] = palette['red']
palette['Sextupole'] = palette['green']
palette['Octupole'] = palette['green']
palette['Multipole'] = palette['green']

# Variables required for the old aperture.py
STYLE_BEND_HATCH = '/'
STYLE_QUAD_HATCH = ''


class ManzoniMatplotlibArtist(_MatplotlibArtist):
    """A matplotlib implementation of a `ZgoubiPlot` artist."""

    def __init__(self,
                 ax=None,
                 with_boxes: bool = True,
                 with_centers: bool = False,
                 tracks_color: str = 'b',
                 **kwargs):
        """
        Args:
            param ax: the matplotlib ax used for plotting. If None it will be created with `init_plot` (kwargs are
            forwarded).
            with_boxes: draw the body of each elements
            with_frames: draw the entry and exit frames of each elements
            with_centers: draw the center of each polar coordinate elements
            tracks_color: color for the plotting of tracks
            kwargs: forwarded to `MatplotlibPlot` and to `init_plot`.
        """
        super().__init__(ax=ax, **kwargs)
        self._with_centers = with_centers
        self._tracks_color = tracks_color

        if ax is None:
            self.init_plot(**kwargs)
        else:
            self._ax = ax
        self._ax2 = self._ax.twinx()
        self._ax2.set_ylim([0, 1])
        self._ax2.axis('off')

    @property
    def tracks_color(self):
        """
        The color for the rendering of the tracks.

        Returns:
            color as a string
        """
        return self._tracks_color

    @tracks_color.setter
    def tracks_color(self, color: str):
        self._tracks_color = color

    @property
    def ax(self):
        """Current Matplotlib ax.

        Returns:
            the Matplotlib ax.
        """
        return self._ax

    @property
    def ax2(self):
        """

        Returns:

        """
        return self._ax2

    @property
    def figure(self):
        """Current Matplotlib figure.

        Returns:
            the Matplotlib figure.
        """
        return self._fig

    @ax.setter
    def ax(self, ax):
        self._ax = ax

    def init_plot(self, figsize=(12, 8), subplots=111):
        """
        Initialize the Matplotlib figure and ax.

        Args:
            subplots: number of subplots
            figsize: figure size
        """
        self._fig = plt.figure(figsize=figsize)
        self._ax = self._fig.add_subplot(subplots)

    def plot(self, *args, **kwargs):
        """Proxy for matplotlib.pyplot.plot

        Same as `matplotlib.pyplot.plot`, forwards all arguments.
        """
        self._ax.plot(*args, **kwargs)

    # THIS IS THE OLD common.py
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

    @staticmethod
    def beamline_get_ticks_locations(o):
        return list(o['AT_CENTER'].values)

    @staticmethod
    def beamline_get_ticks_labels(o):
        return list(o.index)

    @staticmethod
    def draw_beamline(ax, bl):

        offset = 1.15
        ax2 = ax.twinx()
        ax2.set_yticks([])
        ax2.set_ylim([0, 1])
        ax2.hlines(offset, 0, bl['L'].sum().m_as('m'), clip_on=False, colors='black', lw=1)
        for i, e in bl.query("CLASS=='Sbend' or CLASS=='Rbend'").iterrows():
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
                patches.Rectangle(
                    (e['AT_ENTRY'], offset - 0.05 + focusing * 0.02),
                    e['L'].m_as('m'),
                    .1,
                    hatch='',
                    facecolor='blue',
                    clip_on=False,
                )
            )
        for i, e in bl.query("CLASS=='Sextupole' or CLASS=='Quadrupole' or CLASS=='Multipole'").iterrows():
            fc = 'g'
            ax2.add_patch(
                patches.Rectangle(
                    (e['AT_ENTRY'], offset - 0.05),
                    e['L'].m_as('m'),
                    .1,
                    hatch='',
                    facecolor=palette[e['CLASS']],
                    ec=palette[e['CLASS']],
                    clip_on=False,
                )
            )
        for i, e in bl.query("CLASS=='RectangularCollimator'").iterrows():
            fc = 'g'
            ax2.add_patch(
                patches.Rectangle(
                    (e['AT_ENTRY'], offset - 0.05),
                    e['L'].m_as('m'),
                    .1,
                    hatch='',
                    facecolor='gold',
                    ec='gold',
                    clip_on=False,
                )
            )
        for i, e in bl.query("CLASS=='Degrader'").iterrows():
            fc = 'g'
            ax2.add_patch(
                patches.Rectangle(
                    (e['AT_ENTRY'], offset - 0.05),
                    e['L'].m_as('m'),
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
                patches.Rectangle(
                    (e['AT_ENTRY'], offset - 0.05),
                    e['L'].m_as('m'),
                    .1,
                    hatch='',
                    facecolor='grey',
                    ec='grey',
                    clip_on=False,
                )
            )

    def prepare(self, ax, bl, **kwargs):

        bl_short = bl.reset_index()
        bl_short = bl_short[[not a for a in bl_short['NAME'].str.contains("DRIFT")]]
        bl_short = bl_short.set_index("NAME")

        ticks_locations_short = self.beamline_get_ticks_locations(bl_short)
        ticks_labels_short = self.beamline_get_ticks_labels(bl_short)
        ticks_locations = self.beamline_get_ticks_locations(bl)
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
        ax.grid(True, alpha=0.25)

        if kwargs.get('print_label', True):
            ax2 = ax.twiny()
            ax2.set_xlim([ticks_locations[0], ticks_locations[-1]])
            ax2.get_xaxis().set_tick_params(direction='out')
            ax2.tick_params(axis='both', which='major')
            ax2.tick_params(axis='x', labelsize=8)
            plt.setp(ax2.xaxis.get_majorticklabels(), rotation=-90)
            ax2.xaxis.set_major_locator(FixedLocator(ticks_locations_short))
            ax2.xaxis.set_major_formatter(FixedFormatter(ticks_labels_short))

        if kwargs.get("size_arrows", False):
            ax.set_yticklabels([str(abs(x)) for x in ax.get_yticks()])
            ax.annotate('', xy=(-0.103, 0.97), xytext=(-0.103, 0.75),
                        arrowprops=dict(arrowstyle="->", color='k'), xycoords=ax.transAxes)
            ax.annotate('', xy=(-0.103, 0.25), xycoords='axes fraction', xytext=(-0.103, 0.03),
                        arrowprops=dict(arrowstyle="<-", color='k'))
            ax.text(-0.126, 0.86, "Vertical", fontsize=7, rotation=90, transform=ax.transAxes)
            ax.text(-0.126, 0.22, "Horizontal", fontsize=7, rotation=90, transform=ax.transAxes)

        if kwargs.get("with_beamline", False):
            self.draw_beamline(ax, bl)

    @staticmethod
    def filled_plot(ax, x, y0, y, c, fill=False, **kwargs):
        ax.plot(x, y, '.', markersize=0,
                markerfacecolor=c, markeredgecolor=c, color=c, **kwargs)
        if fill:
            ax.fill_between(x, y0, y, facecolor=c, linewidth=0.0, edgecolor=c, **kwargs)

    # THIS IS THE OLD statistics.py
    class StatisticException(Exception):
        """Exception raised for errors in the beam plotting module."""

        def __init__(self, m):
            self.message = m

    @staticmethod
    def histogram_fit(data, bounds_binning=50, verbose=False, model=GaussianModel):
        """ All models are available on https://lmfit.github.io/lmfit-py/builtin_models.html#lmfit.models"""
        y, bin_edges = _np.histogram(data, density=False, bins=bounds_binning)
        x = (bin_edges[:-1] + bin_edges[1:]) / 2
        result = model().fit(data=y, x=x, center=data.mean(), sigma=data.std())
        if verbose:
            print(result.fit_report())
        return x, result

    # THIS IS THE OLD beam.py
    class BeamPlottingException(Exception):
        """Exception raised for errors in the beam plotting module."""

        def __init__(self, m):
            self.message = m

    @staticmethod
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

    @staticmethod
    def rotation_angle(eval, evec):
        if eval[0] * evec[0, 0] > eval[1] * evec[0, 1]:
            return _np.degrees(_np.arctan(evec[1, 0] / evec[0, 0]))
        else:
            return _np.degrees(_np.arctan(evec[1, 1] / evec[0, 1]))

    def phase_space(self, fig, data, **kwargs):
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
        X_mean = _np.mean(x)
        Y_mean = _np.mean(y)

        xlims = [X_mean - 1.2 * _np.max(_np.abs(x)), X_mean + 1.2 * _np.max(_np.abs(x))]
        ylims = [Y_mean - 1.2 * _np.max(_np.abs(y)), Y_mean + 1.2 * _np.max(_np.abs(y))]

        CV = _np.cov(x, y)
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
        h, xedges, yedges = _np.histogram2d(y, x, bins=(ybins, xbins))
        # Plot the beam data
        _ = (ax_global.imshow(h, extent=[xmin, xmax, ymin, ymax],
                              interpolation='nearest', origin='lower',
                              aspect='auto', cmap=kwargs.get('cmap', 'gist_gray_r')))
        # Plot the beam plot contours
        contourcolor = 'black'

        ang = self.rotation_angle(Eval, Evec)
        ang_annotx = _np.cos(_np.radians(ang))
        ang_annoty = _np.sin(_np.radians(ang))
        if kwargs.get('draw_ellipse', True):
            ax_global.add_patch(self.ellipse(X_r[0], Y_r[0], ang, X_mean, Y_mean))
            ax_global.add_patch(self.ellipse(X_r[1], Y_r[1], ang, X_mean, Y_mean))
            ax_global.add_patch(self.ellipse(X_r[2], Y_r[2], ang, X_mean, Y_mean))
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
        xbins = _np.arange(xmin, xmax, (xmax - xmin) / nxbins)
        ybins = _np.arange(ymin, ymax, (ymax - ymin) / nybins)
        # Plot the histograms
        nx, binsx, _ = ax_histx.hist(x, bins=xbins, color='blue', histtype='step')
        ny, binsy, _ = ax_histy.hist(y, bins=ybins, orientation='horizontal', color='red', histtype='step')
        # Set up the histogram limits
        ax_global.set_xlim(xlims)
        ax_global.set_ylim(ylims)
        ax_histx.set_xlim(xlims)
        ax_histy.set_ylim(ylims)
        if kwargs.get('gaussian_fit', False):
            limx = _np.arange(xlims[0], xlims[1], (xlims[1] - xlims[0]) / nxbins)
            bin_centerx, fitresults_x = self.histogram_fit(x, bounds_binning=limx,
                                                           verbose=kwargs.get("verbose", False))
            ax_histx.plot(bin_centerx, fitresults_x.best_fit, 'k--', linewidth=1)
            limy = _np.arange(ylims[0], ylims[1], (ylims[1] - ylims[0]) / nybins)
            bin_centery, fitresults_y = self.histogram_fit(y, bounds_binning=limy,
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

    def phase_space_d(self, ax_global, ax_histx, ax_histy, ax_tab, beam_o_df, elt, dim):
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

        dico_plane = {'X': 0, 'PX': 1, 'Y': 2, 'PY': 3}

        x = 1e3 * _np.array(beam_o_df['BEAM_OUT'][elt][:, dico_plane[dim[0]]])
        y = 1e3 * _np.array(beam_o_df['BEAM_OUT'][elt][:, dico_plane[dim[1]]])

        X_mean = _np.mean(x)
        Y_mean = _np.mean(y)

        xlims = [X_mean - _np.max(_np.abs(x)), X_mean + _np.max(_np.abs(x))]
        ylims = [Y_mean - _np.max(_np.abs(y)), Y_mean + _np.max(_np.abs(y))]

        CV = _np.cov(x, y)
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
        h, xedges, yedges = _np.histogram2d(y, x, bins=(ybins, xbins))
        # Plot the beam data
        _ = (ax_global.imshow(h, extent=[xmin, xmax, ymin, ymax],
                              interpolation='nearest', origin='lower',
                              aspect='auto', cmap='gist_gray_r'))
        # Plot the beam plot contours
        ang = self.rotation_angle(Eval, Evec)
        ang_annotx = _np.cos(_np.radians(ang));
        ang_annoty = _np.sin(_np.radians(ang))
        ax_global.add_patch(self.ellipse(X_r[0], Y_r[0], ang, X_mean, Y_mean, color='red', label='$1\\sigma$'))
        ax_global.add_patch(self.ellipse(X_r[1], Y_r[1], ang, X_mean, Y_mean, color='blue', label='$2\\sigma$'))
        ax_global.add_patch(self.ellipse(X_r[2], Y_r[2], ang, X_mean, Y_mean, color='green', label='$3\\sigma$'))
        ax_global.legend()

        # Plot the axes labels
        ax_global.set_xlabel(dim[0] + ' ' + unit_col_0, fontsize=18)
        ax_global.set_ylabel(dim[1] + ' ' + unit_col_1, fontsize=18)
        # Set up the histogram bins
        xbins = _np.arange(xmin, xmax, (xmax - xmin) / nxbins)
        ybins = _np.arange(ymin, ymax, (ymax - ymin) / nybins)
        # Plot the histograms
        nx, binsx, _ = ax_histx.hist(x, bins=xbins, color='blue', histtype='step')
        ny, binsy, _ = ax_histy.hist(y, bins=ybins, orientation='horizontal', color='red', histtype='step')
        # Set up the histogram limits
        ax_global.set_xlim(xlims)
        ax_global.set_ylim(ylims)
        ax_histx.set_xlim(xlims)
        ax_histy.set_ylim(ylims)

        limx = _np.arange(xlims[0], xlims[1], (xlims[1] - xlims[0]) / nxbins)
        bin_centerx, fitresults_x = self.histogram_fit(x, bounds_binning=limx, verbose=False)
        ax_histx.plot(bin_centerx, fitresults_x.best_fit, 'k--', linewidth=1)
        limy = _np.arange(ylims[0], ylims[1], (ylims[1] - ylims[0]) / nybins)
        bin_centery, fitresults_y = self.histogram_fit(y, bounds_binning=limy, verbose=False)
        ax_histy.plot(fitresults_y.best_fit, bin_centery, 'k--', linewidth=1)

        mean_x = str(round(x.mean(), 3))
        std_x = str(round(x.std(), 3))
        median_x = str(round(_np.median(x), 3))
        mean_y = str(round(y.mean(), 3))
        std_y = str(round(y.std(), 3))
        median_y = str(round(_np.median(y), 3))

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

    def five_spot_map(self, bl_track0, bl_track1, bl_track2, bl_track3, bl_track4, figsize=(8, 8)):
        fig = plt.figure(figsize=figsize)
        try:
            # Order the 5 simulations outputs
            center = _pd.DataFrame(bl_track0.line['BEAM']['ISO'].distribution['X'])
            blcorner = _pd.DataFrame(bl_track1.line['BEAM']['ISO'].distribution['X'])
            brcorner = _pd.DataFrame(bl_track2.line['BEAM']['ISO'].distribution['X'])
            tlcorner = _pd.DataFrame(bl_track3.line['BEAM']['ISO'].distribution['X'])
            trcorner = _pd.DataFrame(bl_track4.line['BEAM']['ISO'].distribution['X'])
            datax = _pd.concat([center, blcorner, brcorner, tlcorner, trcorner], ignore_index=True)

            center = _pd.DataFrame(bl_track0.line['BEAM']['ISO'].distribution['Y'])
            blcorner = _pd.DataFrame(bl_track1.line['BEAM']['ISO'].distribution['Y'])
            brcorner = _pd.DataFrame(bl_track2.line['BEAM']['ISO'].distribution['Y'])
            tlcorner = _pd.DataFrame(bl_track3.line['BEAM']['ISO'].distribution['Y'])
            trcorner = _pd.DataFrame(bl_track4.line['BEAM']['ISO'].distribution['Y'])
            datay = _pd.concat([center, blcorner, brcorner, tlcorner, trcorner], ignore_index=True)
        except:
            print('Error, one simulation is missing. You have to put the 5 simulations outputs in the functions')

        # Create the 2D histogram with the 5 spots.
        _ = plt.hist2d(datax['X'].values, datay['Y'].values, bins=400, cmap='gist_gray_r')
        data_X = []
        data_X.append(1e3 * bl_track0.line['BEAM']['ISO'].std['X'])
        data_X.append(1e3 * bl_track1.line['BEAM']['ISO'].std['X'])
        data_X.append(1e3 * bl_track2.line['BEAM']['ISO'].std['X'])
        data_X.append(1e3 * bl_track3.line['BEAM']['ISO'].std['X'])
        data_X.append(1e3 * bl_track4.line['BEAM']['ISO'].std['X'])

        data_S = []
        data_S.append(100 * _np.abs(bl_track0.line['BEAM']['ISO'].std['X'] - bl_track0.line['BEAM']['ISO'].std['Y']) / (
                bl_track0.line['BEAM']['ISO'].std['X'] + bl_track0.line['BEAM']['ISO'].std['Y']))
        data_S.append(100 * _np.abs(bl_track1.line['BEAM']['ISO'].std['X'] - bl_track1.line['BEAM']['ISO'].std['Y']) / (
                bl_track1.line['BEAM']['ISO'].std['X'] + bl_track1.line['BEAM']['ISO'].std['Y']))
        data_S.append(100 * _np.abs(bl_track2.line['BEAM']['ISO'].std['X'] - bl_track2.line['BEAM']['ISO'].std['Y']) / (
                bl_track2.line['BEAM']['ISO'].std['X'] + bl_track2.line['BEAM']['ISO'].std['Y']))
        data_S.append(100 * _np.abs(bl_track3.line['BEAM']['ISO'].std['X'] - bl_track3.line['BEAM']['ISO'].std['Y']) / (
                bl_track3.line['BEAM']['ISO'].std['X'] + bl_track3.line['BEAM']['ISO'].std['Y']))
        data_S.append(100 * _np.abs(bl_track4.line['BEAM']['ISO'].std['X'] - bl_track4.line['BEAM']['ISO'].std['Y']) / (
                bl_track4.line['BEAM']['ISO'].std['X'] + bl_track4.line['BEAM']['ISO'].std['Y']))

        data_T = []
        data_T.append(_np.degrees(
            self.ellipse_angle_of_rotation(self.fitEllipse(1e3 * bl_track0.line['BEAM']['ISO'].distribution['X'],
                                                           1e3 * bl_track0.line['BEAM']['ISO'].distribution[
                                                               'Y']))))
        data_T.append(_np.degrees(
            self.ellipse_angle_of_rotation(self.fitEllipse(1e3 * bl_track1.line['BEAM']['ISO'].distribution['X'],
                                                           1e3 * bl_track1.line['BEAM']['ISO'].distribution[
                                                               'Y']))))
        data_T.append(_np.degrees(
            self.ellipse_angle_of_rotation(self.fitEllipse(1e3 * bl_track2.line['BEAM']['ISO'].distribution['X'],
                                                           1e3 * bl_track2.line['BEAM']['ISO'].distribution[
                                                               'Y']))))
        data_T.append(_np.degrees(
            self.ellipse_angle_of_rotation(self.fitEllipse(1e3 * bl_track3.line['BEAM']['ISO'].distribution['X'],
                                                           1e3 * bl_track3.line['BEAM']['ISO'].distribution[
                                                               'Y']))))
        data_T.append(_np.degrees(
            self.ellipse_angle_of_rotation(self.fitEllipse(1e3 * bl_track4.line['BEAM']['ISO'].distribution['X'],
                                                           1e3 * bl_track4.line['BEAM']['ISO'].distribution[
                                                               'Y']))))

        plt.annotate(rf"$\sigma_X = {round(data_X[1], 2)}$mm \n $S = {round(data_S[1], 2)}$% \n $\phi={round(data_T[1], 2)} °$"
                     , xy=(- 0.1 + 0.03, 0.1), xytext=(- 0.1 + 0.03, 0.1),
                     horizontalalignment='center', verticalalignment='center',
                     fontsize=12)
        plt.annotate(rf"$\sigma_X = {round(data_X[4], 2)}$mm \n $S = {round(data_S[4], 2)}$% \n $\phi={round(data_T[4], 2)} °$"
                     , xy=(0.1 - 0.03, 0.1), xytext=(0.1 - 0.03, 0.1),
                     horizontalalignment='center', verticalalignment='center',
                     fontsize=12)
        plt.annotate(rf"$\sigma_X = {round(data_X[2], 2)}$mm \n $S = {round(data_S[2], 2)}$% \n $\phi={round(data_T[2], 2)} °$"
                     , xy=(- 0.1 + 0.03, -0.1), xytext=(- 0.1 + 0.03, -0.1),
                     horizontalalignment='center', verticalalignment='center',
                     fontsize=12)
        plt.annotate(rf"$\sigma_X = {round(data_X[3], 2)}$mm \n $S = {round(data_S[3], 2)}$% \n $\phi={round(data_T[3], 2)} °$"
                     , xy=(0.1 - 0.03, -0.1), xytext=(0.1 - 0.03, -0.1),
                     horizontalalignment='center', verticalalignment='center',
                     fontsize=12)
        plt.annotate(rf"$\sigma_X = {round(data_X[0], 2)}$mm \n $S = {round(data_S[0], 2)}$% \n $\phi={round(data_T[0], 2)} °$"
                     , xy=(0.0 + 0.03, 0.0), xytext=(0.0 + 0.03, 0.0),
                     horizontalalignment='center', verticalalignment='center',
                     fontsize=12)

    @staticmethod
    def halo(distribution, dimensions=['X', 'Y', 'PX', 'PY']):
        """Return a dataframe containing the 1st, 5th, 95th and 99th percentiles of each dimensions."""

        halo = _pd.concat([
            distribution[dimensions].quantile(0.01),
            distribution[dimensions].quantile(0.05),
            distribution[dimensions].quantile(1.0 - 0.842701),
            distribution[dimensions].quantile(0.842701),
            distribution[dimensions].quantile(0.95),
            distribution[dimensions].quantile(0.99)
        ], axis=1).rename(columns={0.01: '1%',
                                   0.05: '5%',
                                   1.0 - 0.842701: '20%',
                                   0.842701: '80%',
                                   0.95: '95%',
                                   0.99: '99%'
                                   }
                          )
        return halo

    @staticmethod
    def compute_halo(beam_o_df, element, percentile, dimensions=['X', 'Y', 'PX', 'PY']):
        """Return a dataframe containing the 1st, 5th, 95th and 99th percentiles of each dimensions."""

        data = _pd.DataFrame(columns=['X', 'Y', 'PX', 'PY'])
        data['X'] = beam_o_df['BEAM_OUT'][element][:, 0]
        data['Y'] = beam_o_df['BEAM_OUT'][element][:, 2]
        data['PX'] = beam_o_df['BEAM_OUT'][element][:, 1]
        data['PY'] = beam_o_df['BEAM_OUT'][element][:, 3]

        halo = data[dimensions].quantile(percentile)

        return halo

    # THIS IS THE OLD aperture.py
    @staticmethod
    def xy_from_string(a, i, c):
        def convert(x):
            try:
                converted = float(x)
            except ValueError:
                try:
                    converted = float(c.get(x))
                except TypeError:
                    converted = 1.0
            return converted

        if len(str(a).strip('[]').strip('{}').split(',')) >= int(i) + 1:
            return convert(str(a).strip('[]').strip('{}').split(',')[int(i)])
        elif len(str(a).strip('[]').strip('{}').split(',')) > 0:
            return convert(str(a).strip('[]').strip('{}').split(',')[0])
        else:
            return _np.inf

    @staticmethod
    def draw_chamber(ax, e):
        ax.add_patch(
            patches.Rectangle(
                (e['AT_ENTRY'], (e['APERTURE_UP'])),  # (x,y)
                e['L'].m_as('m'),  # width
                1000 * e['CHAMBER_UP'],  # height
                hatch='', facecolor=palette['base01']
            )
        )
        ax.add_patch(
            patches.Rectangle(
                (e['AT_ENTRY'], -e['APERTURE_DOWN']),  # (x,y)
                e['L'].m_as('m'),  # width
                -1000 * e['CHAMBER_UP'],  # height
                hatch='', facecolor=palette['base01']
            )
        )

    def draw_quad(self, ax, e):
        ax.add_patch(
            patches.Rectangle(
                (e['AT_ENTRY'], e['APERTURE_UP'] + e['CHAMBER_UP']),  # (x,y)
                e['L'].m_as('m'),  # width
                100,  # height
                hatch=STYLE_QUAD_HATCH, facecolor=palette['quad']
            )
        )

        ax.add_patch(
            patches.Rectangle(
                (e['AT_ENTRY'], -e['APERTURE_DOWN'] - e['CHAMBER_UP']),  # (x,y)
                e['L'].m_as('m'),  # width
                -100,  # height
                hatch=STYLE_QUAD_HATCH, facecolor=palette['quad']
            )
        )
        self.draw_chamber(ax, e)

    @staticmethod
    def draw_coll(ax, e, plane):
        # if 'PIPE' not in e:
        #    return
        # if not _np.isnan(e['PIPE']):
        #    return
        ax.add_patch(
            patches.Rectangle(
                (e['AT_ENTRY'], e['APERTURE_UP']),  # (x,y)
                e['L'].m_as('m'),  # width
                100,  # height
                facecolor=palette['coll']
            )
        )

        ax.add_patch(
            patches.Rectangle(
                (e['AT_ENTRY'], -e['APERTURE_DOWN']),  # (x,y)
                e['L'].m_as('m'),  # width
                -100,  # height
                facecolor=palette['coll']
            )
        )

    def draw_bend(self, ax, e):
        tmp = e['APERTURE_UP'] + e['CHAMBER_UP']
        ax.add_patch(
            patches.Rectangle(
                (e['AT_ENTRY'], tmp if tmp < 55 else 55),  # (x,y)
                e['L'].m_as('m'),  # width
                100,  # height
                hatch=STYLE_BEND_HATCH, facecolor=palette['bend']
            )
        )
        tmp = -e['APERTURE_DOWN'] - e['CHAMBER_UP']
        ax.add_patch(
            patches.Rectangle(
                (e['AT_ENTRY'], tmp if abs(tmp) < 55 else -55),  # (x,y)
                e['L'].m_as('m'),  # width
                -100,  # height
                hatch=STYLE_BEND_HATCH, facecolor=palette['bend']
            )
        )
        self.draw_chamber(ax, e)

    @staticmethod
    def fill_aperture(element, context):
        if element.name + '_APERTURE' in context and element['TYPE'] == 'SLITS':
            element['APERTURE'] = context[element.name + '_APERTURE']
        if element['TYPE'] == 'COLLIMATOR' and \
                element['PLUG'] == 'APERTURE' and \
                element['APERTURE'] is not None and \
                _np.isnan(element['APERTURE']):
            element['APERTURE'] = element['CIRCUIT']
        return element

    def aperture(self, ax, bl, **kwargs):

        if 'APERTURE' not in bl:
            return

        planes = kwargs.get('plane', 'both')

        if planes == 'X':
            index = 0
        elif planes == 'Y':
            index = 1

        bl['APERTURE_UP'] = bl['APERTURE'].apply(
            lambda a: a[index].m_as('mm')
        )
        bl['APERTURE_DOWN'] = bl['APERTURE'].apply(
            lambda a: a[index].m_as('mm')
        )

        if 'CHAMBER' not in bl:
            bl['CHAMBER'] = 0

        bl['CHAMBER_UP'] = bl['CHAMBER'].apply(
            lambda a: a
        )
        bl['CHAMBER_DOWN'] = bl['CHAMBER'].apply(
            lambda a: a
        )

        bl.query("CLASS == 'Quadrupole'").apply(lambda e: self.draw_quad(ax, e), axis=1)
        bl.query("CLASS == 'Sbend'").apply(lambda e: self.draw_bend(ax, e), axis=1)
        bl.query("CLASS == 'Rbend'").apply(lambda e: self.draw_bend(ax, e), axis=1)
        bl.query("CLASS == 'RectangularCollimator'").apply(lambda e: self.draw_coll(ax, e, planes), axis=1)
        bl.query("CLASS == 'CircularCollimator'").apply(lambda e: self.draw_coll(ax, e, planes), axis=1)

    # THIS IS THE OLD bpm.py
    @staticmethod
    def draw_bpm_size(ax, s, x):
        ax.add_patch(
            patches.Rectangle(
                (s - 0.05, -x),
                0.1,
                2 * x,
            )
        )

    def bpm(self, ax, bl, **kwargs):
        """TODO."""
        if kwargs.get('plane') is None:
            raise Exception("'plane' argument must be provided.")
        bl.line[bl.line[f"BPM_STD_{kwargs.get('plane')}"].notnull()].apply(
            lambda x: self.draw_bpm_size(ax, x['AT_CENTER'], x[f"BPM_STD_{kwargs.get('plane')}"]),
            axis=1
        )

    # THIS IS THE OLD losses.py
    def losses(self, ax, bl, beam_o_df, **kwargs):
        """Plot the losses from a beamline tracking computation and a context."""
        losses_palette = kwargs.get("palette", palette)
        bl['n_particles'] = list(map(len, beam_o_df['BEAM_OUT']))
        init = bl.query("TYPE == TYPE").drop_duplicates(subset='AT_CENTER', keep='first').iloc[0]['n_particles']
        transmission = bl.query("TYPE == TYPE").drop_duplicates(subset='AT_CENTER', keep='first').apply(
            lambda r: _pd.Series({
                'S': r['AT_EXIT'],
                'T': r['n_particles'] / init
            }), axis=1)
        ax2 = ax.twinx()
        self.prepare(ax2, bl)
        self.prepare(ax, bl, with_beamline=kwargs.get("with_beamline", False))
        ticks_locations = self.beamline_get_ticks_locations(bl)
        ax2.set_ylabel(r'T ($\%$)')
        ax2.yaxis.label.set_color(losses_palette['green'])
        ax2.grid(True)
        if kwargs.get('log', False):
            ax2.semilogy(transmission['S'], 100 * transmission['T'], 's-', color=losses_palette['green'])
            ax2.set_ylim([min(100 * transmission['T']), 100])
        else:
            ax2.yaxis.set_major_locator(MultipleLocator(10))
            ax2.set_ylim([0, 100])
            ax2.plot(transmission['S'], 100 * transmission['T'], 's-', color=losses_palette['green'])
        ax.set_xlim([ticks_locations[0], ticks_locations[-1]])
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.set_ylabel(r'Losses ($\%$)')
        ax.yaxis.label.set_color(losses_palette['magenta'])
        ax.bar(transmission['S'], -100 * transmission['T'].diff(), 0.125, alpha=0.7,
               edgecolor=losses_palette['magenta'],
               color=losses_palette['magenta'],
               error_kw=dict(ecolor=losses_palette['base02'], lw=1, capsize=2, capthick=1))
        ax.set_ylim([0, 100 * transmission['T'].diff().abs().max() + 5.0])

        if kwargs.get('with_current'):
            ax2.plot(bl['AT_EXIT'], bl['CURRENT'], 'k-')

        return transmission, ax, ax2

    # THIS IS THE OLD scattering.py
    @staticmethod
    def draw_slab(ax, e):
        materials_colors = {
            'graphite': 'g',
            'beryllium': 'r',
            'water': 'b',
            'lexan': 'y',
        }
        ax.add_patch(
            patches.Rectangle(
                (e['AT_ENTRY'], -1),  # (x,y)
                e['LENGTH'],  # width
                2,  # height
                hatch='', facecolor=materials_colors.get(e['MATERIAL'], 'k')
            )
        )

    @staticmethod
    def draw_measuring_plane(ax, e):
        ax.add_patch(
            patches.Rectangle(
                (e['AT_ENTRY'] - 0.005, -1),  # (x,y)
                0.01,  # width
                2,  # height
                hatch='', facecolor='k'
            )
        )

    def scattering(self, ax, bl, **kwargs):
        bl.line.query("TYPE == 'slab'").apply(lambda e: self.draw_slab(ax, e), axis=1)
        bl.line.query("TYPE == 'mp'").apply(lambda e: self.draw_measuring_plane(ax, e), axis=1)

    # THIS IS THE OLD summary.py
    def tracking_summary(self, bl=None, fig=None):
        if bl is None:
            raise Exception("`bl` cannot be None.")
        if fig is None:
            fig = plt.figure(figsize=(12, 21))
        ax = fig.add_subplot(311)
        self.prepare(ax, bl, size_arrows=False)
        self.aperture(ax, bl, plane='X')
        self.tracking(ax, bl, plane='X')

        ax = fig.add_subplot(312)
        self.prepare(ax, bl, size_arrows=False)
        self.aperture(ax, bl, plane='Y')
        self.tracking(ax, bl, plane='Y')

        ax = fig.add_subplot(313)
        self.prepare(ax, bl, size_arrows=False)
        self.losses(ax, bl)

        plt.tight_layout()

    def summary(self, bl, beam_o_df, element='DRIFT_ISO', fig=None):
        if fig is None:
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
        self.phase_space_d(ax_global, ax_histx, ax_histy, ax_tab, beam_o_df, element, dim)

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
        self.phase_space_d(ax_global, ax_histx, ax_histy, ax_tab, beam_o_df, element, dim)

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
        self.phase_space_d(ax_global, ax_histx, ax_histy, ax_tab, beam_o_df, element, dim)

        # Define the left top corner block
        n_height = (height + space + hheight) / 3
        n_width = width + space + hwidth
        bottom_h = n_bottom - 1.5 * space + n_height + space
        bottom_hh = bottom_h + n_height + space
        # Set up the geometry of the three plots
        rect_beam_x = [left, bottom_hh, n_width, n_height]
        rect_beam_y = [left, bottom_h, n_width, n_height]
        rect_trans = [left, n_bottom - 1.5 * space, n_width, n_height]
        # Make the three plots
        ax_beam_x = fig.add_axes(rect_beam_x)  # beam plot
        ax_beam_y = fig.add_axes(rect_beam_y)  # x histogram
        ax_trans = fig.add_axes(rect_trans)  # y histogram

        self.prepare(ax_beam_x, bl)
        self.aperture(ax_beam_x, bl, plane='X')
        self.tracking(ax_beam_x, bl, beam_o_df, plane='X', halo=True, halo99=True, std=True, mean=True)

        self.prepare(ax_beam_y, bl, print_label=False)
        self.aperture(ax_beam_y, bl, plane='Y')
        self.tracking(ax_beam_y, bl, beam_o_df, plane='Y', halo=True, halo99=True, std=True, mean=True)

        self.prepare(ax_trans, bl, print_label=False)
        losses_ = self.losses(ax_trans, bl, beam_o_df, log=False)

        ax_beam_x.set_ylabel("Horizontal beam size [mm]")
        ax_beam_y.set_ylabel("Vertical beam size [mm]")
        ax_beam_x.set_xticklabels([])
        ax_beam_y.set_xticklabels([])
        ax_beam_x.set_xlabel('')
        ax_beam_y.set_xlabel('')
        ax_beam_x.set_ylim([-40, 40])
        ax_beam_y.set_ylim([-40, 40])

        return fig

    # THIS SI THE OLD survey.py
    @staticmethod
    def survey_iba(ax, bl, **kwargs):
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Y (m)")
        ax.set_xlim([0, _np.max(bl.line['X']) / 1000])
        ax.set_ylim([0, _np.max(bl.line['Y']) / 1000])
        for index, row in bl.line.iterrows():
            if _pd.notnull(row['X']) and _pd.notnull(row['Y']):
                if kwargs.get("labels", False):
                    ax.annotate(index,
                                xy=(row['X'] / 1000, row['Y'] / 1000),
                                xytext=(row['X'] / 1000 + 0.5, row['Y'] / 1000 + 0.5),
                                arrowprops=dict(arrowstyle="->", facecolor='black', shrinkA=50, shrinkB=5000, ),
                                size=3,
                                horizontalalignment='left',
                                verticalalignment='bottom',
                                clip_on=True
                                )
                if row['CLASS'] == 'QUADRUPOLE':
                    ax.add_patch(
                        patches.Rectangle(
                            (row['X'] / 1000, row['Y'] / 1000),  # (x,y)
                            0.20,  # width
                            0.20,  # height
                            facecolor='#268bd2',
                            edgecolor='#268bd2'
                        )
                    )
                elif row['CLASS'] == 'RBEND' or row['CLASS'] == 'SBEND':
                    ax.add_patch(
                        patches.Rectangle(
                            (row['X'] / 1000, row['Y'] / 1000),  # (x,y)
                            0.250,  # width
                            0.250,  # height
                            facecolor='#dc322f',
                            edgecolor='#dc322f'
                        )
                    )
                else:
                    ax.add_patch(
                        patches.Rectangle(
                            (row['X'] / 1000, row['Y'] / 1000),  # (x,y)
                            0.250,  # width
                            0.250,  # height
                            facecolor='#657b83',
                            edgecolor='#657b83'
                        )
                    )
        plt.axes().set_aspect('equal', 'datalim')

    ELEMENT_COLORS = {
        'QUADRUPOLE': '#cb4b16',  # Solarized 'orange'
        'SEXTUPOLE': 'g',
        'OCTUPOLE': 'g',
        'MULTIPOLE': 'g'
    }

    def survey_madx(self, ax, bl):
        x_edges = [
            _np.min(list(map(_np.abs, [bl.line['Z'].max(), bl.line['Z'].min()]))),
            _np.max(list(map(_np.abs, [bl.line['Z'].max(), bl.line['Z'].min()])))
        ]
        y_edges = [
            _np.min(list(map(_np.abs, [bl.line['X'].max(), bl.line['X'].min()]))),
            _np.max(list(map(_np.abs, [bl.line['X'].max(), bl.line['X'].min()])))
        ]
        edges = [
            _np.max([x_edges[0], y_edges[0]]),
            _np.max([x_edges[1], y_edges[1]])
        ]
        padding = _np.max(edges) * 0.1
        ax.set_xlim([edges[0] - padding, edges[1] + padding])
        ax.set_ylim([edges[0] - padding, edges[1] + padding])
        tmp = -90
        for index, row in bl.line.iterrows():
            if row['KEYWORD'] == 'SBEND' or row['KEYWORD'] == 'MARKER':
                if row['KEYWORD'] == 'SBEND':
                    r = row['L'] / row['ANGLE']
                    if row.get('MAIN_BEND'):
                        c = 'r' if row['MAIN_BEND'] is True else 'b'
                    else:
                        if row['ANGLE'] > 0:
                            c = 'r'
                        elif row['ANGLE'] < 0:
                            c = 'b'
                        else:
                            c = 'gray'
                if row['KEYWORD'] == 'MARKER':
                    r = 0.1
                    c = 'b'
                theta = -row['THETA']
                angle = -row['ANGLE'] / _np.pi * 180.0
                centre = [(row['Z'] - _np.sin(theta) * r), -(row['X'] - _np.cos(theta) * r)]
                theta1 = min(tmp, tmp - angle)
                theta2 = max(tmp, tmp - angle)
                tmp = theta1 if angle > 0 else theta2
                w = patches.Wedge(centre, r + 0.2, theta1, theta2, width=0.4, alpha=1, facecolor=c, ec=c)
            if row['KEYWORD'] == 'DRIFT':
                w = patches.Rectangle(
                    (
                        row['Z'] - row['L'] * _np.cos(row['THETA']) - 0.2 * _np.sin(row['THETA']),
                        -row['X'] + row['L'] * _np.sin(row['THETA']) - 0.2 * _np.cos(row['THETA'])
                    ),
                    row['L'],
                    0.4,
                    angle=_np.degrees(-row['THETA']),
                    alpha=0.2,
                    facecolor='y',
                    ec='y',
                    hatch=''
                )
            if row['KEYWORD'] in ('MULTIPOLE', 'QUADRUPOLE', 'SEXTUPOLE', 'OCTUPOLE'):
                w = patches.Rectangle(
                    (
                        row['Z'] - row['L'] * _np.cos(row['THETA']) - 0.2 * _np.sin(row['THETA']),
                        -row['X'] + row['L'] * _np.sin(row['THETA']) - 0.2 * _np.cos(row['THETA'])
                    ),
                    row['L'],
                    0.4,
                    angle=_np.degrees(-row['THETA']),
                    alpha=1.0,
                    facecolor=self.ELEMENT_COLORS[row['KEYWORD']],
                    ec=self.ELEMENT_COLORS[row['KEYWORD']],
                    hatch=''
                )
            ax.add_patch(w)

    def survey(self, ax, bl, style='madx', **kwargs):
        if style == 'iba_survey':
            self.survey_iba(ax, bl, **kwargs)
        elif style == 'madx':
            self.survey_madx(ax, bl)
        else:
            raise Exception("Style not supported.")

    # THIS IS THE OLD tracking.py
    def tracking(self, ax, bl, beam_o_df, mean=False, std=False, halo=True, **kwargs):
        """Plot the beam envelopes from tracking data."""
        if kwargs.get("plane") is None:
            raise Exception("Plane (plane='X' or plane='Y') must be specified.")

        plane = kwargs.get("plane")
        tracking_palette = kwargs.get("palette", palette)
        if plane is None:
            raise Exception("The 'plane' keyword argument must be set to 'X' or 'Y'.")
        halo_99 = kwargs.get("halo_99")
        std_bpm = kwargs.get("std_bpm", False)

        dico_plane = {'X': 0, 'PX': 1, 'Y': 2, 'PY': 3}

        bl = bl.reset_index()

        t = bl.apply(lambda r: _pd.Series({
            'NAME': r['NAME'],
            'S': r[kwargs.get("reference_plane", 'AT_EXIT')],
            '1%': 1000 * (
                        self.compute_halo(beam_o_df, r['NAME'], 0.023, plane) - self.compute_halo(beam_o_df, r['NAME'],
                                                                                                  0.5,
                                                                                                  plane)) if halo_99 else 0.0,
            '5%': 1000 * (
                        self.compute_halo(beam_o_df, r['NAME'], 0.159, plane) - self.compute_halo(beam_o_df, r['NAME'],
                                                                                                  0.5,
                                                                                                  plane)) if halo else 0.0,
            '95%': 1000 * (
                        self.compute_halo(beam_o_df, r['NAME'], 0.841, plane) - self.compute_halo(beam_o_df, r['NAME'],
                                                                                                  0.5,
                                                                                                  plane)) if halo else 0.0,
            '99%': 1000 * (
                        self.compute_halo(beam_o_df, r['NAME'], 0.977, plane) - self.compute_halo(beam_o_df, r['NAME'],
                                                                                                  0.5,
                                                                                                  plane)) if halo_99 else 0.0,
            'mean': 1000 * beam_o_df['BEAM_OUT'][r['NAME']][:, dico_plane[plane]].mean() if mean else 0.0,
            'std': 1000 * beam_o_df['BEAM_OUT'][r['NAME']][:, dico_plane[plane]].std() if std else 0.0,
        }), axis=1)

        t.set_index('NAME')

        if t['S'].count == 0:
            return

        if halo:
            self.filled_plot(ax, t['S'], t['5%'], t['95%'], tracking_palette[plane], True, alpha=0.3)
            self.filled_plot(ax, t['S'], t['5%'], t['95%'], tracking_palette[plane], True, alpha=0.3)
            if halo_99:
                self.filled_plot(ax, t['S'], t['1%'], t['99%'], tracking_palette[plane], True, alpha=0.3)

        if std:
            ax.plot(t['S'], t['mean'] + t['std'],
                    '^-',
                    color=tracking_palette[plane],
                    markeredgecolor=tracking_palette[plane],
                    markersize=2,
                    linewidth=1
                    )
            ax.plot(t['S'], t['mean'] - t['std'],
                    'v-',
                    color=tracking_palette[plane],
                    markeredgecolor=tracking_palette[plane],
                    markersize=2,
                    linewidth=1
                    )

        if std_bpm:
            # Adjustment to avoid plotting zero values where no BPM is present
            t.loc[t.std_bpm == 0, 'std_bpm'] = -1000
            ax.errorbar(t['S'] - 0.05, t['std_bpm'], xerr=0.1, yerr=t['std_bpm_err'],
                        fmt='none',
                        elinewidth=2.0,
                        linewidth=0.0,
                        color=tracking_palette['green'])
            ax.errorbar(t['S'] - 0.05, -t['std_bpm'], xerr=0.1, yerr=t['std_bpm_err'],
                        fmt='none',
                        elinewidth=2.0,
                        linewidth=0.0,
                        color=tracking_palette['green'])

        if mean:
            ax.plot(t['S'], t['mean'],
                    '*-',
                    color=tracking_palette[plane],
                    markeredgecolor=tracking_palette[plane],
                    markersize=2,
                    linewidth=1,
                    label=kwargs.get("label")
                    )

    # THIS IS THE OLD twiss.py
    def twiss(self, ax, bl, **kwargs):
        """Plot the Twiss beam envelopes from a beamline Twiss computation and a context."""
        context = kwargs.get('context', {})
        bl = bl.line

        bl['XMAXMONO'] = _np.sqrt(bl['BETX'] * context['EMITX'])
        bl['XMAX'] = _np.sqrt(bl['XMAXMONO'] ** 2 + (context['DPP'] * bl['DX']) ** 2)
        bl['YMAX'] = _np.sqrt(bl['BETY'] * context['EMITY'])
        p = kwargs.get('plane', None)
        cx = kwargs.get('color', 'X')
        cy = kwargs.get('color', 'Y')
        if p is None:
            self.filled_plot(ax, bl['S'], 0, -1000 * bl['XMAXMONO'], palette[cx], True, alpha=0.8)
            self.filled_plot(ax, bl['S'], 0, -1000 * bl['XMAX'], palette[cx], True, alpha=0.4)
            self.filled_plot(ax, bl['S'], 0, -2 * 1000 * bl['XMAX'], palette[cx], True, alpha=0.2)
            self.filled_plot(ax, bl['S'], 0, 1000 * bl['YMAX'], palette[cy], True, alpha=0.4)
            self.filled_plot(ax, bl['S'], 0, 2 * 1000 * bl['YMAX'], palette[cy], True, alpha=0.2)
        elif p == 'X':
            self.filled_plot(ax, bl['S'], 0, -1000 * bl['XMAXMONO'], palette[cx], True, alpha=0.8)
            self.filled_plot(ax, bl['S'], 0, -1000 * bl['XMAX'], palette[cx], True, alpha=0.4)
            self.filled_plot(ax, bl['S'], 0, -2 * 1000 * bl['XMAX'], palette[cx], True, alpha=0.2)
            self.filled_plot(ax, bl['S'], 0, 1000 * bl['XMAXMONO'], palette[cx], True, alpha=0.8)
            self.filled_plot(ax, bl['S'], 0, 1000 * bl['XMAX'], palette[cx], True, alpha=0.4)
            self.filled_plot(ax, bl['S'], 0, 2 * 1000 * bl['XMAX'], palette[cx], True, alpha=0.2)
        else:
            self.filled_plot(ax, bl['S'], 0, 1000 * bl['YMAX'], palette[cy], True, alpha=0.4)
            self.filled_plot(ax, bl['S'], 0, 2 * 1000 * bl['YMAX'], palette[cy], True, alpha=0.2)
            self.filled_plot(ax, bl['S'], 0, -1000 * bl['YMAX'], palette[cy], True, alpha=0.4)
            self.filled_plot(ax, bl['S'], 0, -2 * 1000 * bl['YMAX'], palette[cy], True, alpha=0.2)

    def beta(self, ax, bl, **kwargs):
        """Plot the Twiss beta functions."""
        if kwargs.get('ptc', False):
            self.twiss_function_plot(ax, bl, ['BETA'], x='11', y='22', **kwargs)
        else:
            self.twiss_function_plot(ax, bl, ['BET'], **kwargs)

    def alpha(self, ax, bl, **kwargs):
        """Plot the Twiss alpha functions."""
        self.twiss_function_plot(ax, bl, ['ALF'], **kwargs)

    def dispersion(self, ax, bl, planes='both', rel_beta=1, **kwargs):
        """Plot the dispersion functions."""
        if kwargs.get('ptc', True):
            self.twiss_function_plot(ax, bl, ['DISP'], ptc=True, planes=planes)
        else:
            # Caution: MAD-X dispersion is affected by relativistic factors
            # See section 1.7.4 of the MAD-X user guide
            ax.plot(bl.line['S'], rel_beta * bl.line['DX'], color=palette['X'])
            ax.plot(bl.line['S'], rel_beta * bl.line['DY'], color=palette['Y'])

    def phase_advance(self, ax, bl, **kwargs):
        """Plot the phase advance."""
        self.twiss_function_plot(ax, bl, ['MU'], kwargs.get('ptc', False))

    @staticmethod
    def twiss_function_plot(ax, bl, functions, planes='both', **kwargs):
        bl = bl.line

        if kwargs.get('ptc', False):
            x = kwargs.get('x', '1')
            y = kwargs.get('y', '3')
        else:
            x = 'X'
            y = 'Y'
        for f in functions:
            if planes == 'both' or planes == 'X':
                ax.plot(bl['S'], bl[f + x], color=palette['X'])
            if planes == 'both' or planes == 'Y':
                ax.plot(bl['S'], bl[f + y], color=palette['Y'])
