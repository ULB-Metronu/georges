"""Plotly plotting module for Manzoni.

TODO
"""

from __future__ import annotations

import logging

import pandas as _pd
import numpy as _np
from ..manzoni.observers import Observer as _Observer
from ..manzoni.observers import MeanObserver as _MeanObserver
from ..manzoni.observers import SigmaObserver as _SigmaObserver
from ..manzoni.observers import BeamObserver as _BeamObserver
from ..manzoni.observers import LossesObserver as _LossesObserver
from ..manzoni.observers import TwissObserver as _TwissObserver
from georges_core.vis import PlotlyArtist as _PlotlyArtist
from georges_core.vis.artist import PALETTE

palette = PALETTE['solarized']
palette['both'] = palette['base03']
palette['X'] = palette['cyan']
palette['Y'] = palette['orange']
palette['XP'] = palette['red']
palette['YP'] = palette['green']


class BeamPlottingException(Exception):
    """Exception raised for errors in the beam plotting module."""

    def __init__(self, m):
        self.message = m


class ManzoniPlotlyArtist(_PlotlyArtist):
    """A matplotlib implementation of a `Matplotlib` artist."""

    def __init__(self,
                 tracks_color: str = 'b',
                 **kwargs):
        """
        Args:
            param ax: the matplotlib ax used for plotting. If None it will be created with `init_plot` (kwargs are
            forwarded).
            tracks_color: color for the plotting of tracks
            kwargs: forwarded to `ZgoubiPlot` and to `init_plot`.
        """
        super().__init__(**kwargs)
        self._tracks_color = tracks_color

    def tracking(self,
                 observer: _Observer = None,
                 plane: str = 'X',
                 fill_between: bool = False,
                 mean: bool = True,
                 std: bool = False,
                 halo: bool = True,
                 **kwargs):
        """
        Plot the beam envelopes from tracking data.
        Args:
            observer: Observer used for the tracking
            plane:
            fill_between:
            mean:
            std:
            halo:
            **kwargs:

        Returns:

        """
        tracking_palette = kwargs.get("palette", palette)
        df_observer = observer.to_df()

        if isinstance(observer, _MeanObserver):
            self.scatter(x=_np.hstack([0, df_observer['AT_EXIT'].apply(lambda e: e.m_as('m')).values]),
                         y=_np.hstack([df_observer.iloc[0][f"BEAM_IN_{plane}"],
                                       df_observer.iloc[:][f"BEAM_OUT_{plane}"].values]) * 1000,
                         marker={'symbol': 4, 'color': tracking_palette[plane], 'size': 7},
                         line={'width': 1},
                         mode='markers+lines',
                         showlegend=False,
                         name=kwargs.get('name', plane)
                         )
            self.layout['xaxis']['title'] = "S (m)"
            self.layout['yaxis']['title'] = "Mean position (mm)"
        elif isinstance(observer, _SigmaObserver):
            x = _np.hstack([0, df_observer['AT_EXIT'].apply(lambda e: e.m_as('m')).values])

            if plane == 'both':
                y0 = _np.hstack([df_observer.iloc[0][f"BEAM_IN_X"],
                                 df_observer.iloc[:][f"BEAM_OUT_X"].values]) * 1000
                y1 = _np.hstack([df_observer.iloc[0][f"BEAM_IN_Y"],
                                 df_observer.iloc[:][f"BEAM_OUT_Y"].values]) * 1000
            else:
                y0 = _np.hstack([df_observer.iloc[0][f"BEAM_IN_{plane}"],
                                 df_observer.iloc[:][f"BEAM_OUT_{plane}"].values]) * 1000
                y1 = y0
            y = [y1, -y0]

            self.scatter(x=x,
                         y=y[0],
                         marker={'symbol': 4, 'color': tracking_palette[plane], 'size': 7},
                         line={'width': 1, 'color': tracking_palette[plane]},
                         mode='markers+lines',
                         name=kwargs.get('name', plane),
                         showlegend=False,
                         )

            self.scatter(x=x,
                         y=y[1],
                         marker={'symbol': 4, 'color': tracking_palette[plane], 'size': 7},
                         line={'width': 1, 'color': tracking_palette[plane]},
                         mode='markers+lines',
                         showlegend=False,
                         name=kwargs.get('name', plane),
                         fill='tonexty' if fill_between else None,
                         fillcolor=tracking_palette[plane]
                         )

            self.layout['xaxis']['title'] = "S (m)"
            self.layout['yaxis']['title'] = "Beam Size (mm)"

            if plane == 'both':
                self.fig['data'][0]['line']['color'] = tracking_palette['Y']
                self.fig['data'][0]['marker']['color'] = tracking_palette['Y']
                self.fig['data'][1]['line']['color'] = tracking_palette['X']
                self.fig['data'][1]['marker']['color'] = tracking_palette['X']
                ticks_val = _np.arange(_np.min(10 * _np.floor(y[1] / 10)), _np.max(10 * _np.ceil(y[0] / 10)) + 10, 10)
                self.layout['yaxis']['tickvals'] = ticks_val
                self.layout['yaxis']['ticktext'] = [_np.abs(d) for d in ticks_val]

                # Optimized for a general layout
                self.layout['yaxis']['title'] = f"<--------------   Beam size (mm)   -------------->"
                self.layout['annotations'] = [{'x': -0.075,
                                               'y': 0.65 + self.layout['height'] / 6000,
                                               'xref': 'paper',
                                               'yref': 'paper',
                                               'text': "Vertical",
                                               'textangle': -90},
                                              {'x': -0.075,
                                               'y': self.layout['height'] / 6000,
                                               'xref': 'paper',
                                               'yref': 'paper',
                                               'text': "Horizontal",
                                               'textangle': -90}]

        elif isinstance(observer, _BeamObserver):

            dico_plane = {'X': 0, 'PX': 1, 'Y': 2, 'PY': 3}
            t = df_observer.apply(lambda r: _pd.Series({
                'S': r['AT_EXIT'].m_as('m'),
                'mean': 1000 * r['BEAM_OUT'][:, dico_plane[plane]].mean() if mean else 0.0,
                'std': 1000 * r['BEAM_OUT'][:, dico_plane[plane]].std() if std else 0.0,
                '1%': 1000 * (self.compute_halo(r['BEAM_OUT'][:, dico_plane[plane]], 0.023)
                              - self.compute_halo(r['BEAM_OUT'][:, dico_plane[plane]], 0.5)),
                '5%': 1000 * (self.compute_halo(r['BEAM_OUT'][:, dico_plane[plane]], 0.159)
                              - self.compute_halo(r['BEAM_OUT'][:, dico_plane[plane]], 0.5)),
                '95%': 1000 * (self.compute_halo(r['BEAM_OUT'][:, dico_plane[plane]], 0.841)
                               - self.compute_halo(r['BEAM_OUT'][:, dico_plane[plane]], 0.5)),
                '99%': 1000 * (self.compute_halo(r['BEAM_OUT'][:, dico_plane[plane]], 0.977)
                               - self.compute_halo(r['BEAM_OUT'][:, dico_plane[plane]], 0.5)),
            }), axis=1)

            if df_observer.iloc[0]['BEAM_IN'] is not None:
                data_entry = df_observer.iloc[0]
                t0 = _pd.DataFrame(data={"S": [data_entry['AT_ENTRY'].m_as('m')],
                                         'mean': [1000 * data_entry['BEAM_IN'][
                                                         :, dico_plane[plane]].mean() if mean else 0.0],
                                         'std': [
                                             1000 * data_entry['BEAM_IN'][:, dico_plane[plane]].std() if std else 0.0],
                                         '1%': [1000 * (self.compute_halo(data_entry['BEAM_IN'][:, dico_plane[plane]],
                                                                          0.023)
                                                        - self.compute_halo(data_entry['BEAM_IN'][:, dico_plane[plane]],
                                                                            0.5))],
                                         '5%': [1000 * (self.compute_halo(data_entry['BEAM_IN'][:, dico_plane[plane]],
                                                                          0.159)
                                                        - self.compute_halo(data_entry['BEAM_IN'][:, dico_plane[plane]],
                                                                            0.5))],
                                         '95%': [1000 * (self.compute_halo(data_entry['BEAM_IN'][:, dico_plane[plane]],
                                                                           0.841)
                                                         - self.compute_halo(
                                                     data_entry['BEAM_IN'][:, dico_plane[plane]],
                                                     0.5))],
                                         '99%': [1000 * (self.compute_halo(data_entry['BEAM_IN'][:, dico_plane[plane]],
                                                                           0.977)
                                                         - self.compute_halo(
                                                     data_entry['BEAM_IN'][:, dico_plane[plane]],
                                                     0.5))]
                                         },
                                   index=["Start"])
                t = _pd.concat([t0, t])

            if t['S'].count == 0:
                return

            if halo:
                self.filled_plot(self, t['S'], t['mean'] + t['5%'], t['mean'] + t['95%'],
                                 tracking_palette[plane], fill=True, alpha=0.3)
                self.filled_plot(self, t['S'], t['mean'] + t['1%'], t['mean'] + t['99%'],
                                 tracking_palette[plane], fill=True, alpha=0.3)
            #
            if mean:
                self.scatter(x=t['S'],
                             y=t['mean'],
                             marker={'symbol': 2, 'color': tracking_palette[plane], 'size': 7},
                             line={'width': 1},
                             mode='markers+lines',
                             showlegend=False,
                             name='mean',
                             )

            if std:
                self.scatter(x=t['S'],
                             y=t['mean'] + t['std'],
                             marker={'symbol': 2, 'color': tracking_palette[plane], 'size': 7},
                             line={'width': 1},
                             mode='markers+lines',
                             showlegend=False,
                             name='std',
                             )

                self.scatter(x=t['S'],
                             y=t['mean'] - t['std'],
                             marker={'symbol': 2, 'color': tracking_palette[plane], 'size': 7},
                             line={'width': 1},
                             mode='markers+lines',
                             showlegend=False,
                             name='std',
                             )

            self.layout['xaxis']['title'] = "S (m)"
            self.layout['yaxis']['title'] = "Beam Size (mm)"

        elif isinstance(observer, _LossesObserver):
            raise BeamPlottingException(f"Use method vis.ManzoniPlotlyArtist().losses to plot losses.")

        else:
            raise BeamPlottingException(f"No plotting method for {observer} is implemented")

    def losses(self,
               observer: _LossesObserver = None,
               log_scale: bool = False,
               **kwargs):
        """
        Plot the losses along the beamline
        Args:
            observer: Observer used for the tracking
            log_scale: Log scale for transmission
        """

        if not isinstance(observer, _LossesObserver):
            raise BeamPlottingException("The observer must be a LossesObserver.")

        losses_palette = kwargs.get("palette", palette)
        df_observer = observer.to_df()
        exit = df_observer['AT_EXIT'].apply(lambda e: e.m_as('m'))

        self.layout['xaxis']['title'] = "S (m)"
        self.layout['yaxis']['title'] = r"Losses (%)"
        self.layout['yaxis']['titlefont'] = {'color': losses_palette['magenta']}

        self.bar(x=exit,
                 y=df_observer['LOSSES'],
                 marker={'color': losses_palette['magenta']},
                 width=0.125,
                 name='Losses',
                 showlegend=False)

        max_val = (df_observer['LOSSES']).abs().max()
        self.layout['yaxis']['dtick'] = _np.ceil(max_val / 10)
        self.layout['yaxis']['range'] = [0, max_val + 5.0]

        self.add_secondary_axis(title=r'T (%)')
        self.layout['yaxis2']['titlefont'] = {'color': losses_palette['green']}

        # TODO USE Transmission, shift and compute ?
        init = df_observer.iloc[0]['PARTICLES_IN']
        global_transmission = 100 * (df_observer['PARTICLES_OUT'].values / init)

        if log_scale:
            self.scatter(x=_np.hstack([0, exit.values]),
                         y=_np.hstack([100, global_transmission]),
                         mode='lines+markers',
                         marker={'symbol': 17, 'color': losses_palette['green']},
                         yaxis='y2',
                         name='Transmission',
                         showlegend=False)
            self.layout['yaxis2']['type'] = 'log'
            self.layout['yaxis2']['range'] = [_np.floor(_np.log(_np.min(global_transmission))), 2]
            self.layout['yaxis2']['tickvals'] = 10**(_np.arange(_np.floor(_np.log(_np.min(global_transmission))), 3, 1))

        else:
            self.scatter(x=_np.hstack([0, exit.values]),
                         y=_np.hstack([100, global_transmission]),
                         yaxis='y2',
                         mode='lines+markers',
                         marker={'symbol': 1, 'color': losses_palette['green']},
                         name='Transmission',
                         showlegend=False)
            self.layout['yaxis2']['dtick'] = 10
            self.layout['yaxis2']['range'] = [0, 100]

    @staticmethod
    def compute_halo(data, percentile):
        """Return a dataframe containing the 1st, 5th, 95th and 99th percentiles of each dimensions."""
        return _np.quantile(data, percentile)

    @staticmethod
    def filled_plot(ax, x, y0, y, c, fill=False, **kwargs):
        ax.scatter(x=x, y=y, mode='lines', line={'width': 1, 'color': c}, showlegend=False)
        if fill:
            ax.scatter(x=x, y=y0, mode='lines', line={'width': 1, 'color': c}, showlegend=False, fill='tonexty')

    def twiss(self, observer: _TwissObserver = None,):
        logging.error("The method is not yet implemented")

    def phase_space(self, observer: _TwissObserver = None, ):
        logging.error("The method is not yet implemented")
