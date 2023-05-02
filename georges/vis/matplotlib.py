"""Matplotlib plotting module for Manzoni.

TODO
"""

from __future__ import annotations

import logging
from typing import Union

import cpymad.madx
import matplotlib.patches as patches
import matplotlib.ticker as mticker
import numpy as _np
import pandas as _pd
from georges_core.vis import MatplotlibArtist as _MatplotlibArtist
from georges_core.vis.artist import PALETTE
from lmfit.models import GaussianModel

from ..manzoni.observers import BeamObserver as _BeamObserver
from ..manzoni.observers import IbaBpmObserver as _IbaBpmObserver
from ..manzoni.observers import LossesObserver as _LossesObserver
from ..manzoni.observers import MeanObserver as _MeanObserver
from ..manzoni.observers import Observer as _Observer
from ..manzoni.observers import SigmaObserver as _SigmaObserver
from ..manzoni.observers import SymmetryObserver as _SymmetryObserver
from ..manzoni.observers import TwissObserver as _TwissObserver

palette = PALETTE["solarized"]
palette["both"] = palette["base03"]
palette["X"] = palette["cyan"]
palette["Y"] = palette["orange"]
palette["XP"] = palette["red"]
palette["YP"] = palette["green"]
palette["TWISS_X"] = palette["blue"]
palette["TWISS_Y"] = palette["red"]
palette["DX"] = palette["green"]
palette["DY"] = palette["orange"]


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

    def __init__(self, tracks_color: str = "b", **kwargs):
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

    @staticmethod
    def histogram_fit(data, bounds_binning=50, verbose=False, model=GaussianModel):
        """All models are available on https://lmfit.github.io/lmfit-py/builtin_models.html#lmfit.models"""
        y, bin_edges = _np.histogram(data, density=False, bins=bounds_binning)
        x = (bin_edges[:-1] + bin_edges[1:]) / 2
        result = model().fit(data=y, x=x, center=data.mean(), sigma=data.std())
        if verbose:
            print(result.fit_report())
        return x, result

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
        return patches.Ellipse(
            (x0, y0),
            width,
            height,
            angle=angle,
            linewidth=kwargs.get("linewidth", 2),
            fill=kwargs.get("fill", False),
            linestyle=kwargs.get("linestyle", "--"),
            edgecolor=kwargs.get("color", "red"),
            label=kwargs.get("label"),
        )

    @staticmethod
    def rotation_angle(e_val, evec):
        if e_val[0] * evec[0, 0] > e_val[1] * evec[0, 1]:
            return _np.degrees(_np.arctan(evec[1, 0] / evec[0, 0]))
        else:
            return _np.degrees(_np.arctan(evec[1, 1] / evec[0, 1]))

    # Plotting for the tracking
    def tracking(
        self,
        observer: _Observer = None,
        plane: str = "X",
        fill_between: bool = False,
        mean: bool = True,
        std: bool = False,
        halo: bool = True,
        **kwargs,
    ):
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
            self._ax.plot(
                _np.hstack([0, df_observer["AT_EXIT"].apply(lambda e: e.m_as("m")).values]),
                _np.hstack(
                    [
                        df_observer.iloc[0][f"BEAM_IN_{plane}"],
                        df_observer.iloc[:][f"BEAM_OUT_{plane}"].values,
                    ],
                )
                * 1000,
                "^-",
                color=tracking_palette[plane],
                markeredgecolor=tracking_palette[plane],
                markersize=2,
                linewidth=1,
                label=kwargs.get("label", plane),
            )
            self._ax.set_xlabel("S (m)")
            self._ax.set_ylabel("Mean Position (mm)")

        elif isinstance(observer, _SigmaObserver):
            x = _np.hstack([0, df_observer["AT_EXIT"].apply(lambda e: e.m_as("m")).values])
            if plane == "both":
                y0 = (
                    _np.hstack(
                        [
                            df_observer.iloc[0]["BEAM_IN_X"],
                            df_observer.iloc[:]["BEAM_OUT_X"].values,
                        ],
                    )
                    * 1000
                )
                y1 = (
                    _np.hstack(
                        [
                            df_observer.iloc[0]["BEAM_IN_Y"],
                            df_observer.iloc[:]["BEAM_OUT_Y"].values,
                        ],
                    )
                    * 1000
                )
            else:
                y0 = (
                    _np.hstack(
                        [
                            df_observer.iloc[0][f"BEAM_IN_{plane}"],
                            df_observer.iloc[:][f"BEAM_OUT_{plane}"].values,
                        ],
                    )
                    * 1000
                )
                y1 = y0
            y = [y1, -y0]

            if plane == "both":
                label = ""
            else:
                label = kwargs.get("label", plane)
            self._ax.plot(
                x,
                y[0],
                x,
                y[1],
                marker="^",
                color=tracking_palette[plane],
                markeredgecolor=tracking_palette[plane],
                markersize=2,
                linewidth=1,
                label=label,
            )
            if fill_between:
                self._ax.fill_between(
                    x,
                    y[1],
                    y[0],
                    facecolor=tracking_palette[plane],
                    linewidth=0.0,
                    edgecolor=tracking_palette[plane],
                )
            self._ax.set_xlabel("S (m)")
            self._ax.set_ylabel("Beam Size (mm)")

            if plane == "both":
                self._ax.get_lines()[0].set_color(tracking_palette["Y"])
                self._ax.get_lines()[1].set_color(tracking_palette["X"])
                ticks_loc = self._ax.get_yticks().tolist()
                self._ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
                self._ax.set_yticklabels([abs(x) for x in ticks_loc])

                self._ax.annotate(
                    "",
                    xy=(-0.08 / (self._ax.get_figure().get_size_inches()[0] / 10), 0.97),
                    xytext=(-0.08 / (self._ax.get_figure().get_size_inches()[0] / 10), 0.75),
                    arrowprops=dict(arrowstyle="->", color="k"),
                    xycoords=self._ax.transAxes,
                )

                self._ax.annotate(
                    "",
                    xy=(-0.08 / (self._ax.get_figure().get_size_inches()[0] / 10), 0.25),
                    xycoords="axes fraction",
                    xytext=(-0.08 / (self._ax.get_figure().get_size_inches()[0] / 10), 0.03),
                    arrowprops=dict(arrowstyle="<-", color="k"),
                )

                self._ax.text(
                    -0.1 / (self._ax.get_figure().get_size_inches()[0] / 10),
                    0.86,
                    "Vertical",
                    fontsize=9,
                    rotation=90,
                    transform=self._ax.transAxes,
                )
                self._ax.text(
                    -0.1 / (self._ax.get_figure().get_size_inches()[0] / 10),
                    0.11,
                    "Horizontal",
                    fontsize=9,
                    rotation=90,
                    transform=self._ax.transAxes,
                )

        elif isinstance(observer, _BeamObserver):
            dico_plane = {"X": 0, "PX": 1, "Y": 2, "PY": 3}
            t = df_observer.apply(
                lambda r: _pd.Series(
                    {
                        "S": r["AT_EXIT"].m_as("m"),
                        "mean": 1000 * r["BEAM_OUT"][:, dico_plane[plane]].mean() if mean else 0.0,
                        "std": 1000 * r["BEAM_OUT"][:, dico_plane[plane]].std() if std else 0.0,
                        "1%": 1000
                        * (
                            self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 0.023)
                            - self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 0.5)
                        ),
                        "5%": 1000
                        * (
                            self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 0.159)
                            - self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 0.5)
                        ),
                        "20%": 1000
                        * (
                            self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 1 - 0.842701)
                            - self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 0.5)
                        ),
                        "80%": 1000
                        * (
                            self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 0.842701)
                            - self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 0.5)
                        ),
                        "95%": 1000
                        * (
                            self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 0.841)
                            - self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 0.5)
                        ),
                        "99%": 1000
                        * (
                            self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 0.977)
                            - self.compute_halo(r["BEAM_OUT"][:, dico_plane[plane]], 0.5)
                        ),
                    },
                ),
                axis=1,
            )

            if df_observer.iloc[0]["BEAM_IN"] is not None:
                data_entry = df_observer.iloc[0]
                t0 = _pd.DataFrame(
                    data={
                        "S": [data_entry["AT_ENTRY"].m_as("m")],
                        "mean": [
                            1000
                            * data_entry["BEAM_IN"][
                                :,
                                dico_plane[plane],
                            ].mean()
                            if mean
                            else 0.0,
                        ],
                        "std": [
                            1000 * data_entry["BEAM_IN"][:, dico_plane[plane]].std() if std else 0.0,
                        ],
                        "1%": [
                            1000
                            * (
                                self.compute_halo(
                                    data_entry["BEAM_IN"][:, dico_plane[plane]],
                                    0.023,
                                )
                                - self.compute_halo(
                                    data_entry["BEAM_IN"][:, dico_plane[plane]],
                                    0.5,
                                )
                            ),
                        ],
                        "5%": [
                            1000
                            * (
                                self.compute_halo(
                                    data_entry["BEAM_IN"][:, dico_plane[plane]],
                                    0.159,
                                )
                                - self.compute_halo(
                                    data_entry["BEAM_IN"][:, dico_plane[plane]],
                                    0.5,
                                )
                            ),
                        ],
                        "95%": [
                            1000
                            * (
                                self.compute_halo(
                                    data_entry["BEAM_IN"][:, dico_plane[plane]],
                                    0.841,
                                )
                                - self.compute_halo(
                                    data_entry["BEAM_IN"][:, dico_plane[plane]],
                                    0.5,
                                )
                            ),
                        ],
                        "99%": [
                            1000
                            * (
                                self.compute_halo(
                                    data_entry["BEAM_IN"][:, dico_plane[plane]],
                                    0.977,
                                )
                                - self.compute_halo(
                                    data_entry["BEAM_IN"][:, dico_plane[plane]],
                                    0.5,
                                )
                            ),
                        ],
                    },
                    index=["Start"],
                )
                t = _pd.concat([t0, t])

            if t["S"].count == 0:
                return

            if halo:
                self.filled_plot(
                    self._ax,
                    t["S"],
                    t["mean"] + t["5%"],
                    t["mean"] + t["95%"],
                    tracking_palette[plane],
                    fill=True,
                    alpha=0.3,
                )
                self.filled_plot(
                    self._ax,
                    t["S"],
                    t["mean"] + t["1%"],
                    t["mean"] + t["99%"],
                    tracking_palette[plane],
                    fill=True,
                    alpha=0.3,
                )

            if mean:
                self._ax.plot(
                    t["S"],
                    t["mean"],
                    "*-",
                    color=tracking_palette[plane],
                    markeredgecolor=tracking_palette[plane],
                    markersize=2,
                    linewidth=1,
                    label="mean",
                )

            if std:
                self._ax.plot(
                    t["S"],
                    t["mean"] + t["std"],
                    "^-",
                    color=tracking_palette[plane],
                    markeredgecolor=tracking_palette[plane],
                    markersize=2,
                    linewidth=1,
                    label="std",
                )
                self._ax.plot(
                    t["S"],
                    t["mean"] - t["std"],
                    "v-",
                    color=tracking_palette[plane],
                    markeredgecolor=tracking_palette[plane],
                    markersize=2,
                    linewidth=1,
                )

            self._ax.set_xlabel("S (m)")
            self._ax.set_ylabel("Beam Size (mm)")

        elif isinstance(observer, _IbaBpmObserver):
            logging.warning("The plotting method for IbaBpmObserver is not yet implemented.")
            # Adjustment to avoid plotting zero values where no BPM is present
            # if std_bpm:
            #
            #     t.loc[t.std_bpm == 0, 'std_bpm'] = -1000
            #     ax.errorbar(t['S'] - 0.05, t['std_bpm'], xerr=0.1, yerr=t['std_bpm_err'],
            #                 fmt='none',
            #                 elinewidth=2.0,
            #                 linewidth=0.0,
            #                 color=tracking_palette['green'])
            #     ax.errorbar(t['S'] - 0.05, -t['std_bpm'], xerr=0.1, yerr=t['std_bpm_err'],
            #                 fmt='none',
            #                 elinewidth=2.0,
            #                 linewidth=0.0,
            #                 color=tracking_palette['green'])

        elif isinstance(observer, _SymmetryObserver):
            raise BeamPlottingException("Use method vis.ManzoniMatplotlibArtist(ax=ax).symmetry to plot symmetry.")

        elif isinstance(observer, _LossesObserver):
            raise BeamPlottingException("Use method vis.ManzoniMatplotlibArtist(ax=ax).losses to plot losses.")

        else:
            raise BeamPlottingException(f"No plotting method for {observer} is implemented")

    # Plotting for the losses
    def losses(self, observer: _LossesObserver = None, log_scale: bool = False, **kwargs):
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
        data_center = df_observer["AT_CENTER"].apply(lambda e: e.m_as("m"))

        self._ax.set_xlabel(r"S (m)")
        self._ax.set_ylabel(r"Losses ($\%$)")
        self._ax.yaxis.label.set_color(losses_palette["magenta"])
        self._ax.bar(
            data_center,
            df_observer["LOSSES"],
            width=0.125,
            alpha=0.7,
            edgecolor=losses_palette["magenta"],
            color=losses_palette["magenta"],
            align="center",
            error_kw=dict(ecolor=losses_palette["base02"], capsize=2, capthick=1),
        )
        max_val = (df_observer["LOSSES"]).abs().max()
        self._ax.yaxis.set_major_locator(mticker.MultipleLocator(_np.ceil(max_val / 10)))
        self._ax.set_ylim([0, max_val + 5.0])

        ax2 = self._ax.twinx()  # self.ax_2 is reserved for cartouche
        ax2.set_ylabel(r"T ($\%$)")
        ax2.yaxis.label.set_color(losses_palette["green"])
        ax2.grid(True)
        # TODO USE Transmission, shift and compute ?
        init = df_observer.iloc[0]["PARTICLES_IN"]
        global_transmission = 100 * (df_observer["PARTICLES_OUT"].values / init)
        if log_scale:
            ax2.semilogy(
                _np.hstack([0, data_center.values]),
                _np.hstack([100, global_transmission]),
                "s-",
                color=losses_palette["green"],
            )
            ax2.set_ylim([min(global_transmission), 100])
        else:
            ax2.yaxis.set_major_locator(mticker.MultipleLocator(10))
            ax2.set_ylim([0, 100])
            ax2.plot(
                _np.hstack([0, data_center.values]),
                _np.hstack([100, global_transmission]),
                "s-",
                color=losses_palette["green"],
            )

    def symmetry(self, observer: _LossesObserver = None, **kwargs):
        """
        Plot the symmetry of the beam along the beamline

        Args:
            observer: Observer used for the tracking
        """

        if not isinstance(observer, _SymmetryObserver):
            raise BeamPlottingException("The observer must be a SymmetryObserver.")
        symmetry_palette = kwargs.get("palette", palette)

        df_observer = observer.to_df()
        self._ax.plot(
            _np.hstack([0, df_observer["AT_EXIT"].apply(lambda e: e.m_as("m")).values]),
            _np.hstack(
                [
                    df_observer.iloc[0]["SYM_IN"],
                    df_observer.iloc[:]["SYM_OUT"].values,
                ],
            )
            * 100,
            "^-",
            color=symmetry_palette["blue"],
            markeredgecolor=symmetry_palette["blue"],
            markersize=2,
            linewidth=1,
        )
        self._ax.set_xlabel("S (m)")
        self._ax.set_ylabel("Asymmetry ($\%$)")  # noqa: W605

    @staticmethod
    def compute_halo(data, percentile):
        """Return a dataframe containing the 1st, 5th, 95th and 99th percentiles of each dimensions."""
        return _np.quantile(data, percentile)

    @staticmethod
    def filled_plot(ax, x, y0, y, c, fill=False, **kwargs):
        ax.plot(x, y, ".", markersize=0, markerfacecolor=c, markeredgecolor=c, color=c, **kwargs)
        if fill:
            ax.fill_between(x, y0, y, facecolor=c, linewidth=0.0, edgecolor=c, **kwargs)

    # Plotting for the Twiss
    def twiss(
        self,
        observer: _TwissObserver = None,
        with_beta: bool = True,
        with_alpha: bool = False,
        with_dispersion: bool = False,
        tfs_data: Union[_pd.DataFrame, cpymad.madx.Table] = None,
        relativistic_beta: float = 1.0,
        **kwargs,
    ):
        """
        Plot the Twiss function along the beamline

        Args:
            observer: Observer used for the tracking
            with_beta: plot the beta
            with_alpha: plot the alpha
            with_dispersion: plot the dispersion
            tfs_data: if provided, plot the data from MAD-X.
            relativistic_beta (float): Relativistic beta value to scale the dispersion. Default to 1.

        """
        if not isinstance(observer, _TwissObserver):
            raise BeamPlottingException("The observer must be a TwissObserver.")
        if isinstance(tfs_data, cpymad.madx.Table):
            tfs_data = _pd.DataFrame(
                data={
                    "S": tfs_data.s,
                    "BETX": tfs_data.betx,
                    "BETY": tfs_data.bety,
                    "ALFX": tfs_data.alfx,
                    "ALFY": tfs_data.alfy,
                    "DX": tfs_data.dx,
                    "DY": tfs_data.dy,
                },
            )
        df_observer = observer.to_df()
        twiss_palette = kwargs.get("palette", palette)

        if with_beta:
            self._ax.plot(
                _np.hstack([0, df_observer["AT_EXIT"].apply(lambda e: e.m_as("m")).values]),
                _np.hstack(
                    [
                        df_observer.iloc[0]["BETA_IN_X"],
                        df_observer.iloc[:]["BETA_OUT_X"].values,
                    ],
                ),
                "-",
                color=twiss_palette["TWISS_X"],
                linewidth=1,
                label="BETX - Manzoni",
            )

            self._ax.plot(
                _np.hstack([0, df_observer["AT_EXIT"].apply(lambda e: e.m_as("m")).values]),
                _np.hstack(
                    [
                        df_observer.iloc[0]["BETA_IN_Y"],
                        df_observer.iloc[:]["BETA_OUT_Y"].values,
                    ],
                ),
                "-",
                color=twiss_palette["TWISS_Y"],
                linewidth=1,
                label="BETY - Manzoni",
            )
            if tfs_data is not None:
                self._ax.plot(
                    tfs_data["S"].values,
                    tfs_data["BETX"].values,
                    color=twiss_palette["TWISS_X"],
                    markeredgecolor=twiss_palette["TWISS_X"],
                    markersize=4,
                    marker="x",
                    ls="None",
                    label="BETX - MADX",
                )

                self._ax.plot(
                    tfs_data["S"].values,
                    tfs_data["BETY"].values,
                    color=twiss_palette["TWISS_Y"],
                    markeredgecolor=twiss_palette["TWISS_Y"],
                    markersize=4,
                    marker="x",
                    ls="None",
                    label="BETY - MADX",
                )

            self._ax.set_xlabel("S (m)")
            self._ax.set_ylabel(r"$\beta$ (m)")
            max_val = _np.ceil(_np.maximum(df_observer["BETA_OUT_X"].max(), df_observer["BETA_OUT_Y"].max()))
            self._ax.yaxis.set_major_locator(mticker.MultipleLocator(_np.ceil(max_val / 10)))
            self._ax.set_ylim([0, max_val + 5.0])

        if with_alpha:
            self._ax.plot(
                _np.hstack([0, df_observer["AT_EXIT"].apply(lambda e: e.m_as("m")).values]),
                _np.hstack(
                    [
                        df_observer.iloc[0]["ALPHA_IN_X"],
                        df_observer.iloc[:]["ALPHA_OUT_X"].values,
                    ],
                ),
                "-",
                color=twiss_palette["TWISS_X"],
                linewidth=1,
                label="ALPHAX - Manzoni",
            )

            self._ax.plot(
                _np.hstack([0, df_observer["AT_EXIT"].apply(lambda e: e.m_as("m")).values]),
                _np.hstack(
                    [
                        df_observer.iloc[0]["ALPHA_IN_Y"],
                        df_observer.iloc[:]["ALPHA_OUT_Y"].values,
                    ],
                ),
                "-",
                color=twiss_palette["TWISS_Y"],
                linewidth=1,
                label="ALPHAY - Manzoni",
            )
            if tfs_data is not None:
                self._ax.plot(
                    tfs_data["S"].values,
                    tfs_data["ALFX"].values,
                    color=twiss_palette["TWISS_X"],
                    markeredgecolor=twiss_palette["TWISS_X"],
                    markersize=4,
                    marker="x",
                    ls="None",
                    label="ALPHAX - MADX",
                )

                self._ax.plot(
                    tfs_data["S"].values,
                    tfs_data["ALFY"].values,
                    color=twiss_palette["TWISS_Y"],
                    markeredgecolor=twiss_palette["TWISS_Y"],
                    markersize=4,
                    marker="x",
                    ls="None",
                    label="ALPHAY - MADX",
                )

            self._ax.set_xlabel("S (m)")
            self._ax.set_ylabel(r"$\alpha$ (m)")
            max_val = _np.ceil(_np.maximum(df_observer["ALPHA_OUT_X"].max(), df_observer["ALPHA_OUT_Y"].max()))
            min_val = _np.ceil(_np.minimum(df_observer["ALPHA_OUT_X"].min(), df_observer["ALPHA_OUT_Y"].min()))
            self._ax.yaxis.set_major_locator(mticker.MultipleLocator(_np.ceil(max_val / 10)))
            self._ax.set_ylim([min_val - 5.0, max_val + 5.0])

        if with_dispersion:
            self.ax_disp = self._ax.twinx()
            self.ax_disp.plot(
                _np.hstack([0, df_observer["AT_EXIT"].apply(lambda e: e.m_as("m")).values]),
                _np.hstack(
                    [
                        df_observer.iloc[0]["DISP_IN_X"],
                        df_observer.iloc[:]["DISP_OUT_X"].values,
                    ],
                ),
                "-",
                color=twiss_palette["DX"],
                linewidth=1,
                label="DX - Manzoni",
            )

            self.ax_disp.plot(
                _np.hstack([0, df_observer["AT_EXIT"].apply(lambda e: e.m_as("m")).values]),
                _np.hstack(
                    [
                        df_observer.iloc[0]["DISP_IN_Y"],
                        df_observer.iloc[:]["DISP_OUT_Y"].values,
                    ],
                ),
                "-",
                color=twiss_palette["DY"],
                linewidth=1,
                label="DY - Manzoni",
            )
            if tfs_data is not None:
                logging.warning(
                    f"Dispersion from MAD-X is multiplied by the beta relativistic factor {relativistic_beta}.",
                )
                self.ax_disp.plot(
                    tfs_data["S"].values,
                    tfs_data["DX"].values * relativistic_beta,
                    color=twiss_palette["DX"],
                    markeredgecolor=twiss_palette["DX"],
                    markersize=4,
                    marker="x",
                    ls="None",
                    label="DX - MADX",
                )

                self.ax_disp.plot(
                    tfs_data["S"].values,
                    tfs_data["DY"].values * relativistic_beta,
                    color=twiss_palette["DY"],
                    markeredgecolor=twiss_palette["DY"],
                    markersize=4,
                    marker="x",
                    ls="None",
                    label="DY - MADX",
                )

            self.ax_disp.set_xlabel("S (m)")
            self.ax_disp.set_ylabel("Dispersion (m)")
            max_val = _np.ceil(_np.maximum(df_observer["DISP_OUT_X"].max(), df_observer["DISP_OUT_Y"].max()))
            min_val = _np.floor(_np.minimum(df_observer["DISP_OUT_X"].min(), df_observer["DISP_OUT_Y"].min()))
            if max_val > 0:
                self.ax_disp.yaxis.set_major_locator(mticker.MultipleLocator(_np.ceil(max_val / 10)))
            self.ax_disp.set_ylim([min_val - 5.0, max_val + 5.0])
            self.ax_disp.legend()

    # Plotting for the phase space
    def phase_space(
        self,
        observer: _BeamObserver = None,
        element: str = None,
        location: str = "OUT",
        dim=None,
        nbins=None,
        draw_ellipse: bool = True,
    ):
        """

        Args:
            observer:
            element:
            location:
            dim:
            nbins:
            draw_ellipse:
        Returns:

        """

        if nbins is None:
            nbins = [50, 50]

        if dim is None:
            dim = ["X", "Y"]

        if not isinstance(observer, _BeamObserver):
            raise BeamPlottingException("The observer must be a BeamObserver.")

        df_observer = observer.to_df()
        if element is None:
            element = df_observer.iloc[0].name

        data_element = df_observer.loc[element, f"BEAM_{location}"]
        if location == "IN":
            s_position = df_observer.loc[element, "AT_ENTRY"]
        else:
            s_position = df_observer.loc[element, "AT_EXIT"]

        if dim[0] == "X" or dim[0] == "Y":
            unit_col_0 = "[mm]"
        else:
            unit_col_0 = "[mrad]"
        if dim[1] == "X" or dim[1] == "Y":
            unit_col_1 = "[mm]"
        else:
            unit_col_1 = "[mrad]"

        dico_plane = {"X": 0, "PX": 1, "Y": 2, "PY": 3}
        x = 1000 * data_element[:, dico_plane[dim[0]]]
        y = 1000 * data_element[:, dico_plane[dim[1]]]

        self._prepare(sposition=s_position)

        # Main Figure
        xlim = [_np.mean(x) - 5 * _np.std(x), _np.mean(x) + 5 * _np.std(x)]
        ylim = [_np.mean(y) - 5 * _np.std(y), _np.mean(y) + 5 * _np.std(y)]
        h, xedges, yedges = _np.histogram2d(x, y, bins=[nbins[0], nbins[1]], range=[xlim, ylim])
        self._fig.axes[len(self._fig.axes) - 4].imshow(
            h.T,
            extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
            interpolation="nearest",
            origin="lower",
            aspect="auto",
            cmap="gist_gray_r",
        )

        self._fig.axes[len(self._fig.axes) - 4].set_xlabel(f"{dim[0]} {unit_col_0}")
        self._fig.axes[len(self._fig.axes) - 4].set_ylabel(f"{dim[1]} {unit_col_1}")
        if draw_ellipse:
            self.draw_ellipse(x, y)

        self._fig.axes[len(self._fig.axes) - 3].hist(x, bins=nbins[0], color="blue", histtype="step")
        self._fig.axes[len(self._fig.axes) - 2].hist(
            y,
            bins=nbins[1],
            orientation="horizontal",
            color="red",
            histtype="step",
        )
        bin_centerx, fitresults_x = self.histogram_fit(x, bounds_binning=nbins[0], verbose=False)
        self._fig.axes[len(self._fig.axes) - 3].plot(bin_centerx, fitresults_x.best_fit, "k--", linewidth=1)
        bin_centery, fitresults_y = self.histogram_fit(y, bounds_binning=nbins[0], verbose=False)
        self._fig.axes[len(self._fig.axes) - 2].plot(fitresults_y.best_fit, bin_centery, "k--", linewidth=1)

        # Table
        self.draw_table(x, y, dim, unit_col_0, unit_col_1, element)

    def _prepare(self, sposition):
        if len(self._fig.axes) == 1:  # We don't have a cartouche
            space = 0.02
            width, height = 0.25, 0.25
            left, hwidth = 0.05, width / 2
            bottom, hheight = 0.05, height / 2
            bottom_h = left_h = left + width + space

        else:
            bounds = self._ax.get_position().bounds
            space = 0.02
            width, height = 0.8 * (bounds[2] - bounds[0]) - space, 0.8 * (bounds[3] - bounds[1]) - space
            left, hwidth = bounds[0], bounds[2] - width - space
            bottom, hheight = bounds[1], bounds[3] - height - space
            bottom_h = left + height + space
            left_h = left + width + space
            self._ax2.add_patch(
                patches.Rectangle(
                    (sposition.m_as("m"), 1.15 - 0.075),
                    0.01,
                    0.15,
                    hatch="",
                    facecolor="m",
                    clip_on=False,
                ),
            )
            self._ax2.spines["top"].set_visible(False)
            self._ax2.spines["bottom"].set_visible(False)
            self._ax2.spines["right"].set_visible(False)
            self._ax2.spines["left"].set_visible(False)

        self._fig.delaxes(self._ax)
        # Set up the geometry of the three plots
        rect_beam = [left, bottom, width, height]  # dimensions of temp plot
        rect_histx = [left, bottom_h, width, hheight]  # dimensions of x-histogram
        rect_histy = [left_h, bottom, hwidth, height]  # dimensions of y-histogram
        rect_tab = [left_h, bottom_h, hwidth, hheight]  # dimensions of tab
        # Make the three plots
        ax_global = self._fig.add_axes(rect_beam)  # noqa: F841
        ax_histx = self._fig.add_axes(rect_histx)  # x histogram
        ax_histy = self._fig.add_axes(rect_histy)  # y histogram
        ax_histx.set_ylabel("Counts")
        ax_histx.grid(True)
        ax_histy.set_xlabel("Counts")
        ax_histy.grid(True)
        nullfmt = mticker.NullFormatter()
        ax_histx.xaxis.set_major_formatter(nullfmt)
        ax_histy.yaxis.set_major_formatter(nullfmt)
        ax_tab = self._fig.add_axes(rect_tab)  # y histogram
        ax_tab.tick_params(labelbottom="off", labelleft="off", left="off", bottom="off")

    def draw_ellipse(self, x, y):
        [e_val, evec] = _np.linalg.eig(_np.cov(x, y))
        chisq = [2.278868566, 5.991464547, 11.61828598]
        x_r, y_r = [], []
        for i in range(0, len(chisq)):
            x_r.append(2 * (chisq[i] * e_val[0]) ** 0.5)
            y_r.append(2 * (chisq[i] * e_val[1]) ** 0.5)

        ang = self.rotation_angle(e_val, evec)

        self._fig.axes[len(self._fig.axes) - 4].add_patch(
            self.ellipse(
                x_r[0],
                y_r[0],
                ang,
                _np.mean(x),
                _np.mean(y),
                color="red",
                label="$1\\sigma$",
            ),
        )
        self._fig.axes[len(self._fig.axes) - 4].add_patch(
            self.ellipse(
                x_r[1],
                y_r[1],
                ang,
                _np.mean(x),
                _np.mean(y),
                color="blue",
                label="$2\\sigma$",
            ),
        )
        self._fig.axes[len(self._fig.axes) - 4].add_patch(
            self.ellipse(
                x_r[2],
                y_r[2],
                ang,
                _np.mean(x),
                _np.mean(y),
                color="green",
                label="$3\\sigma$",
            ),
        )
        self._fig.axes[len(self._fig.axes) - 4].legend()

    def draw_table(self, x, y, dim, unit_col_0, unit_col_1, element):
        mean_x = f"{_np.round(_np.mean(x), 3)}"
        std_x = f"{_np.round(_np.std(x), 3)}"
        median_x = f"{_np.round(_np.median(x), 3)}"
        mean_y = f"{_np.round(_np.mean(y), 3)}"
        std_y = f"{_np.round(_np.std(y), 3)}"
        median_y = f"{_np.round(_np.median(y), 3)}"

        self._fig.axes[-1].tick_params(labelbottom="off", labelleft="off", left="off", bottom="off")
        self._fig.axes[-1].get_xaxis().set_visible(False)
        self._fig.axes[-1].get_yaxis().set_visible(False)
        x0 = self._fig.axes[-1].get_xlim()[0]
        y0 = self._fig.axes[-1].get_ylim()[0]
        self._fig.axes[-1].axvline(x0 + 0.25, color="k", linewidth=1)
        self._fig.axes[-1].axvline(x0 + 0.5, color="k", linewidth=1)
        self._fig.axes[-1].axvline(x0 + 0.75, color="k", linewidth=1)
        self._fig.axes[-1].axhline(y0 + 0.4, color="k", linewidth=1)
        self._fig.axes[-1].axhline(y0 + 0.8, color="k", linewidth=1)
        self._fig.axes[-1].axhline(y0 + 0.83, color="k", linewidth=1)
        self._fig.axes[-1].annotate(
            element,
            xy=(x0 + 0.125, y0 + 0.9),
            xytext=(x0 + 0.125, y0 + 0.9),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )
        self._fig.axes[-1].annotate(
            "mean",
            xy=(x0 + 0.375, y0 + 0.9),
            xytext=(x0 + 0.375, y0 + 0.9),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )
        self._fig.axes[-1].annotate(
            "std ",
            xy=(x0 + 0.625, y0 + 0.9),
            xytext=(x0 + 0.625, y0 + 0.9),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )
        self._fig.axes[-1].annotate(
            "median",
            xy=(x0 + 0.875, y0 + 0.9),
            xytext=(x0 + 0.875, y0 + 0.9),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )

        xname = f"{dim[0]} {unit_col_0}"
        yname = f"{dim[1]} {unit_col_1}"
        self._fig.axes[-1].annotate(
            xname,
            xy=(x0 + 0.125, 0.6),
            xytext=(x0 + 0.125, 0.6),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=11,
        )
        self._fig.axes[-1].annotate(
            yname,
            xy=(x0 + 0.125, 0.2),
            xytext=(x0 + 0.125, 0.2),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )
        self._fig.axes[-1].annotate(
            mean_x,
            xy=(x0 + 0.375, 0.6),
            xytext=(x0 + 0.375, 0.6),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )
        self._fig.axes[-1].annotate(
            std_x,
            xy=(x0 + 0.625, 0.6),
            xytext=(x0 + 0.625, 0.6),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )
        self._fig.axes[-1].annotate(
            median_x,
            xy=(x0 + 0.875, 0.6),
            xytext=(x0 + 0.875, 0.6),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )
        self._fig.axes[-1].annotate(
            mean_y,
            xy=(x0 + 0.375, 0.2),
            xytext=(x0 + 0.375, 0.2),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )
        self._fig.axes[-1].annotate(
            std_y,
            xy=(x0 + 0.625, 0.2),
            xytext=(x0 + 0.625, 0.2),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )
        self._fig.axes[-1].annotate(
            median_y,
            xy=(x0 + 0.875, 0.2),
            xytext=(x0 + 0.875, 0.2),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )

    # Plotting five spot map
    def five_spot_map(self, bl_track0, bl_track1, bl_track2, bl_track3, bl_track4):
        # TODO
        pass
        # try:
        #     # Order the 5 simulations outputs
        #     center = _pd.DataFrame(bl_track0.line['BEAM']['ISO'].distribution['X'])
        #     blcorner = _pd.DataFrame(bl_track1.line['BEAM']['ISO'].distribution['X'])
        #     brcorner = _pd.DataFrame(bl_track2.line['BEAM']['ISO'].distribution['X'])
        #     tlcorner = _pd.DataFrame(bl_track3.line['BEAM']['ISO'].distribution['X'])
        #     trcorner = _pd.DataFrame(bl_track4.line['BEAM']['ISO'].distribution['X'])
        #     datax = _pd.concat([center, blcorner, brcorner, tlcorner, trcorner], ignore_index=True)
        #
        #     center = _pd.DataFrame(bl_track0.line['BEAM']['ISO'].distribution['Y'])
        #     blcorner = _pd.DataFrame(bl_track1.line['BEAM']['ISO'].distribution['Y'])
        #     brcorner = _pd.DataFrame(bl_track2.line['BEAM']['ISO'].distribution['Y'])
        #     tlcorner = _pd.DataFrame(bl_track3.line['BEAM']['ISO'].distribution['Y'])
        #     trcorner = _pd.DataFrame(bl_track4.line['BEAM']['ISO'].distribution['Y'])
        #     datay = _pd.concat([center, blcorner, brcorner, tlcorner, trcorner], ignore_index=True)
        # except:
        #     print('Error, one simulation is missing. You have to put the 5 simulations outputs in the functions')
        #
        # # Create the 2D histogram with the 5 spots.
        # _ = plt.hist2d(datax['X'].values, datay['Y'].values, bins=400, cmap='gist_gray_r')
        # data_X = []
        # data_X.append(1e3 * bl_track0.line['BEAM']['ISO'].std['X'])
        # data_X.append(1e3 * bl_track1.line['BEAM']['ISO'].std['X'])
        # data_X.append(1e3 * bl_track2.line['BEAM']['ISO'].std['X'])
        # data_X.append(1e3 * bl_track3.line['BEAM']['ISO'].std['X'])
        # data_X.append(1e3 * bl_track4.line['BEAM']['ISO'].std['X'])
        #
        # data_S = []
        # data_S.append(100 * _np.abs(bl_track0.line['BEAM']['ISO'].std['X'] -
        # bl_track0.line['BEAM']['ISO'].std['Y']) / (
        #         bl_track0.line['BEAM']['ISO'].std['X'] + bl_track0.line['BEAM']['ISO'].std['Y']))
        # data_S.append(100 * _np.abs(bl_track1.line['BEAM']['ISO'].std['X'] - bl_track1.line['BEAM']['ISO']
        # .std['Y']) / (
        #         bl_track1.line['BEAM']['ISO'].std['X'] + bl_track1.line['BEAM']['ISO'].std['Y']))
        # data_S.append(100 * _np.abs(bl_track2.line['BEAM']['ISO'].std['X'] - bl_track2.line['BEAM']['ISO']
        # .std['Y']) / (
        #         bl_track2.line['BEAM']['ISO'].std['X'] + bl_track2.line['BEAM']['ISO'].std['Y']))
        # data_S.append(100 * _np.abs(bl_track3.line['BEAM']['ISO'].std['X'] - bl_track3.line['BEAM']['ISO']
        # .std['Y']) / (
        #         bl_track3.line['BEAM']['ISO'].std['X'] + bl_track3.line['BEAM']['ISO'].std['Y']))
        # data_S.append(100 * _np.abs(bl_track4.line['BEAM']['ISO'].std['X'] - bl_track4.line['BEAM']['ISO']
        # .std['Y']) / (
        #         bl_track4.line['BEAM']['ISO'].std['X'] + bl_track4.line['BEAM']['ISO'].std['Y']))
        #
        # data_T = []
        # data_T.append(_np.degrees(
        #     self.ellipse_angle_of_rotation(self.fitEllipse(1e3 * bl_track0.line['BEAM']['ISO'].distribution['X'],
        #                                                    1e3 * bl_track0.line['BEAM']['ISO'].distribution[
        #                                                        'Y']))))
        # data_T.append(_np.degrees(
        #     self.ellipse_angle_of_rotation(self.fitEllipse(1e3 * bl_track1.line['BEAM']['ISO'].distribution['X'],
        #                                                    1e3 * bl_track1.line['BEAM']['ISO'].distribution[
        #                                                        'Y']))))
        # data_T.append(_np.degrees(
        #     self.ellipse_angle_of_rotation(self.fitEllipse(1e3 * bl_track2.line['BEAM']['ISO'].distribution['X'],
        #                                                    1e3 * bl_track2.line['BEAM']['ISO'].distribution[
        #                                                        'Y']))))
        # data_T.append(_np.degrees(
        #     self.ellipse_angle_of_rotation(self.fitEllipse(1e3 * bl_track3.line['BEAM']['ISO'].distribution['X'],
        #                                                    1e3 * bl_track3.line['BEAM']['ISO'].distribution[
        #                                                        'Y']))))
        # data_T.append(_np.degrees(
        #     self.ellipse_angle_of_rotation(self.fitEllipse(1e3 * bl_track4.line['BEAM']['ISO'].distribution['X'],
        #                                                    1e3 * bl_track4.line['BEAM']['ISO'].distribution[
        #                                                        'Y']))))
        #
        # plt.annotate(
        #     rf"$\sigma_X = {round(data_X[1], 2)}$mm \n $S = {round(data_S[1], 2)}$% \n $\phi={round(data_T[1], 2)} °$"
        #     , xy=(- 0.1 + 0.03, 0.1), xytext=(- 0.1 + 0.03, 0.1),
        #     horizontalalignment='center', verticalalignment='center',
        #     fontsize=12)
        # plt.annotate(
        #     rf"$\sigma_X = {round(data_X[4], 2)}$mm \n $S = {round(data_S[4], 2)}$% \n $\phi={round(data_T[4], 2)} °$"
        #     , xy=(0.1 - 0.03, 0.1), xytext=(0.1 - 0.03, 0.1),
        #     horizontalalignment='center', verticalalignment='center',
        #     fontsize=12)
        # plt.annotate(
        #     rf"$\sigma_X = {round(data_X[2], 2)}$mm \n $S = {round(data_S[2], 2)}$% \n $\phi={round(data_T[2], 2)} °$"
        #     , xy=(- 0.1 + 0.03, -0.1), xytext=(- 0.1 + 0.03, -0.1),
        #     horizontalalignment='center', verticalalignment='center',
        #     fontsize=12)
        # plt.annotate(
        #     rf"$\sigma_X = {round(data_X[3], 2)}$mm \n $S = {round(data_S[3], 2)}$% \n $\phi={round(data_T[3], 2)} °$"
        #     , xy=(0.1 - 0.03, -0.1), xytext=(0.1 - 0.03, -0.1),
        #     horizontalalignment='center', verticalalignment='center',
        #     fontsize=12)
        # plt.annotate(
        #     rf"$\sigma_X = {round(data_X[0], 2)}$mm \n $S = {round(data_S[0], 2)}$% \n $\phi={round(data_T[0], 2)} °$"
        #     , xy=(0.0 + 0.03, 0.0), xytext=(0.0 + 0.03, 0.0),
        #     horizontalalignment='center', verticalalignment='center',
        #     fontsize=12)
