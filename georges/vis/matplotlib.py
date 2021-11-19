"""Matplotlib plotting module for Manzoni.

TODO
"""

from __future__ import annotations
import pandas as _pd
import numpy as _np
from ..manzoni.observers import Observer as _Observer
from ..manzoni.observers import MeanObserver as _MeanObserver
from ..manzoni.observers import SigmaObserver as _SigmaObserver
from ..manzoni.observers import BeamObserver as _BeamObserver
# from numpy.linalg import eig
# import matplotlib.pyplot as plt
# import matplotlib.patches as patches
# from matplotlib.ticker import MultipleLocator, NullFormatter
# from lmfit.models import GaussianModel
from georges_core.vis import MatplotlibArtist as _MatplotlibArtist
from georges_core.vis.artist import PALETTE

palette = PALETTE['solarized']
palette['X'] = palette['cyan']
palette['Y'] = palette['orange']
palette['XP'] = palette['red']
palette['YP'] = palette['green']


# THIS IS THE OLD beam.py
class BeamPlottingException(Exception):
    """Exception raised for errors in the beam plotting module."""

    def __init__(self, m):
        self.message = m


class StatisticException(Exception):
    """Exception raised for errors in the beam plotting module."""

    def __init__(self, m):
        self.message = m


class ManzoniMatplotlibArtist(_MatplotlibArtist):
    """A matplotlib implementation of a `Matplotlib` artist."""

    def __init__(self,
                 tracks_color: str = 'b',
                 **kwargs):
        """
        Args:
            param ax: the matplotlib ax used for plotting. If None it will be created with `init_plot` (kwargs are
            forwarded).
            tracks_color: color for the plotting of tracks
            kwargs: forwarded to `MatplotlibPlot` and to `init_plot`.
        """
        super().__init__(**kwargs)
        self._tracks_color = tracks_color

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

    def tracking(self,
                 beamline: _pd.DataFrame = None,
                 observers: _Observer = None,
                 plane: str = 'X',
                 mean: bool = True,
                 std: bool = False,
                 halo: bool = True,
                 **kwargs):
        """
        Plot the beam envelopes from tracking data.
        Args:
            beamline: dataframe of the beamline
            observers: Observer used for the tracking
            plane:
            mean:
            std:
            halo:
            **kwargs:

        Returns:

        """
        # Merge the observers with the beamline
        bl_tracking = _pd.merge(beamline, observers.to_df(), on='NAME')
        tracking_palette = kwargs.get("palette", palette)

        if isinstance(observers, _MeanObserver) or isinstance(observers, _SigmaObserver):
            entry = bl_tracking['AT_ENTRY'].apply(lambda e: e.m_as('m'))
            exit = bl_tracking['AT_EXIT'].apply(lambda e: e.m_as('m'))

            self._ax.plot(entry.values,
                          bl_tracking[f"BEAM_IN_{plane}"]*1000,
                          '^-',
                          color=tracking_palette[plane],
                          markeredgecolor=tracking_palette[plane],
                          markersize=2,
                          linewidth=1
                          )

            self._ax.plot(exit.values,
                          bl_tracking[f"BEAM_OUT_{plane}"]*1000,
                          '^-',
                          color=tracking_palette[plane],
                          markeredgecolor=tracking_palette[plane],
                          markersize=2,
                          linewidth=1
                          )

            self._ax.set_xlabel("S (m)")
            if isinstance(observers, _MeanObserver):
                self._ax.set_ylabel("Mean Position (mm)")
            else:
                self._ax.set_ylabel("Beam Size (mm)")


    def tracking2(self, ax, bl, beam_o_df, mean=False, std=False, halo=True, **kwargs):

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

        # if kwargs.get("size_arrows", False):
        #     ax.set_yticklabels([str(abs(x)) for x in ax.get_yticks()])
        #     ax.annotate('', xy=(-0.103, 0.97), xytext=(-0.103, 0.75),
        #                 arrowprops=dict(arrowstyle="->", color='k'), xycoords=ax.transAxes)
        #     ax.annotate('', xy=(-0.103, 0.25), xycoords='axes fraction', xytext=(-0.103, 0.03),
        #                 arrowprops=dict(arrowstyle="<-", color='k'))
        #     ax.text(-0.126, 0.86, "Vertical", fontsize=7, rotation=90, transform=ax.transAxes)
        #     ax.text(-0.126, 0.22, "Horizontal", fontsize=7, rotation=90, transform=ax.transAxes)

    # THIS IS THE OLD twiss.py
    # TODO do the same thing as in Zgoubidoo

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

# @staticmethod
# def filled_plot(ax, x, y0, y, c, fill=False, **kwargs):
#     ax.plot(x, y, '.', markersize=0,
#             markerfacecolor=c, markeredgecolor=c, color=c, **kwargs)
#     if fill:
#         ax.fill_between(x, y0, y, facecolor=c, linewidth=0.0, edgecolor=c, **kwargs)
#
# @staticmethod
# def halo(distribution, dimensions=['X', 'Y', 'PX', 'PY']):
#     """Return a dataframe containing the 1st, 5th, 95th and 99th percentiles of each dimensions."""
#
#     halo = _pd.concat([
#         distribution[dimensions].quantile(0.01),
#         distribution[dimensions].quantile(0.05),
#         distribution[dimensions].quantile(1.0 - 0.842701),
#         distribution[dimensions].quantile(0.842701),
#         distribution[dimensions].quantile(0.95),
#         distribution[dimensions].quantile(0.99)
#     ], axis=1).rename(columns={0.01: '1%',
#                                0.05: '5%',
#                                1.0 - 0.842701: '20%',
#                                0.842701: '80%',
#                                0.95: '95%',
#                                0.99: '99%'
#                                }
#                       )
#     return halo
#
# @staticmethod
# def compute_halo(beam_o_df, element, percentile, dimensions=['X', 'Y', 'PX', 'PY']):
#     """Return a dataframe containing the 1st, 5th, 95th and 99th percentiles of each dimensions."""
#
#     data = _pd.DataFrame(columns=['X', 'Y', 'PX', 'PY'])
#     data['X'] = beam_o_df['BEAM_OUT'][element][:, 0]
#     data['Y'] = beam_o_df['BEAM_OUT'][element][:, 2]
#     data['PX'] = beam_o_df['BEAM_OUT'][element][:, 1]
#     data['PY'] = beam_o_df['BEAM_OUT'][element][:, 3]
#
#     halo = data[dimensions].quantile(percentile)
#
#     return halo

# # THIS IS THE OLD aperture.py
# @staticmethod
# def xy_from_string(a, i, c):
#     def convert(x):
#         try:
#             converted = float(x)
#         except ValueError:
#             try:
#                 converted = float(c.get(x))
#             except TypeError:
#                 converted = 1.0
#         return converted
#
#     if len(str(a).strip('[]').strip('{}').split(',')) >= int(i) + 1:
#         return convert(str(a).strip('[]').strip('{}').split(',')[int(i)])
#     elif len(str(a).strip('[]').strip('{}').split(',')) > 0:
#         return convert(str(a).strip('[]').strip('{}').split(',')[0])
#     else:
#         return _np.inf
#
# @staticmethod
# def fill_aperture(element, context):
#     if element.name + '_APERTURE' in context and element['TYPE'] == 'SLITS':
#         element['APERTURE'] = context[element.name + '_APERTURE']
#     if element['TYPE'] == 'COLLIMATOR' and \
#             element['PLUG'] == 'APERTURE' and \
#             element['APERTURE'] is not None and \
#             _np.isnan(element['APERTURE']):
#         element['APERTURE'] = element['CIRCUIT']
#     return element
