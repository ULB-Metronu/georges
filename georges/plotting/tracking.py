from .common import filled_plot
from .common import palette as common_palette
import pandas as pd
import numpy as np


def tracking(ax, bl, **kwargs):
    """Plot the beam envelopes from tracking data."""

    plane = kwargs.get("plane")
    palette = kwargs.get("palette", 'X')

    if plane is None:
        raise Exception("The 'plane' keyword argument must be set to 'X' or 'Y'.")
    halo = kwargs.get("halo", True)
    halo_99 = kwargs.get("halo99", halo)
    std = kwargs.get("std", True)
    std_bpm = kwargs.get("std_bpm", False)
    mean = kwargs.get("mean", True)

    t = bl.line.query("BEAM == BEAM").apply(lambda r: pd.Series({
        'S': r[kwargs.get("reference_plane", 'AT_CENTER')],
        '1%': 1000 * r['BEAM'].halo['1%'][plane] if halo_99 else 0.0,
        '5%': 1000 * r['BEAM'].halo['5%'][plane] if halo else 0.0,
        '95%': 1000 * r['BEAM'].halo['95%'][plane] if halo else 0.0,
        '99%': 1000 * r['BEAM'].halo['99%'][plane] if halo_99 else 0.0,
        'mean': 1000 * r['BEAM'].mean[plane] if mean else 0.0,
        'std': 1000 * r['BEAM'].std[plane] if std else 0.0,
        'std_bpm': 1000 * r['BEAM'].std_bpm[plane][0] * int(pd.notnull(r['BPM'])) if 'BPM' in bl.line.columns else 0.0,
        'std_bpm_err': np.max([1.0, 1000 * r['BEAM'].std_bpm[plane][1] * int(pd.notnull(r['BPM'])) if 'BPM' in bl.line.columns else 0.0]),
    }), axis=1)

    if t['S'].count == 0:
        return

    if halo:
        if halo_99:
            filled_plot(ax, t['S'], t['1%'], t['99%'], common_palette[palette], True, alpha=0.3)
        filled_plot(ax, t['S'], t['5%'], t['95%'], common_palette[palette], True, alpha=0.3)
        filled_plot(ax, t['S'], t['mean'] - t['std'], t['mean'] + t['std'], common_palette[palette], True, alpha=0.3)

    if std:
        ax.plot(t['S'], t['mean'] + t['std'], '^-', color=common_palette[palette],
                markeredgecolor=common_palette[palette], markersize=2, linewidth=1)
        ax.plot(t['S'], t['mean'] - t['std'], 'v-', color=common_palette[palette],
                markeredgecolor=common_palette[palette], markersize=2, linewidth=1)

    if std_bpm:
        # Adjustment to avoid plotting zero values where no BPM is present
        t.loc[t.std_bpm == 0, 'std_bpm'] = -1000
        ax.errorbar(t['S'] - 0.05, t['std_bpm'], xerr=0.1, yerr=t['std_bpm_err'],
                    fmt='none',
                    elinewidth=2.0,
                    linewidth=0.0,
                    color=palette['green'])
        ax.errorbar(t['S'] - 0.05, -t['std_bpm'], xerr=0.1, yerr=t['std_bpm_err'],
                    fmt='none',
                    elinewidth=2.0,
                    linewidth=0.0,
                    color=palette['green'])

    if mean:
        ax.plot(t['S'], t['mean'], '*-',
                color=common_palette[plane],
                markeredgecolor=common_palette[plane],
                markersize=2,
                linewidth=1,
                label=kwargs.get("label")
                )


def compare_beamline(ax, bl1, bl2, **kwargs):
    """Compare results obtained with different softwares"""
    if kwargs.get("plane") is None:
        raise Exception("Plane (plane='X' or plane='Y') must be specified.")

    plane = kwargs.get("plane")
    palette = kwargs.get("palette", 'X')

    if plane is None:
        raise Exception("The 'plane' keyword argument must be set to 'X' or 'Y'.")

    std = kwargs.get("std", True)
    mean = kwargs.get("mean", False)

    t1 = bl1.line.query("BEAM == BEAM").apply(lambda r: pd.Series({
        'S': r[kwargs.get("reference_plane_beamline1", 'AT_CENTER')],
        'mean_bl1': 1000 * r['BEAM'].mean[plane],
        'std_bl1': 1000 * r['BEAM'].std[plane]
    }), axis=1)

    t2 = bl2.line.query("BEAM == BEAM").apply(lambda r: pd.Series({
        'S': r[kwargs.get("reference_plane_beamline2", 'AT_CENTER')],
        'mean_bl2': 1000 * r['BEAM'].mean[plane],
        'std_bl2': 1000 * r['BEAM'].std[plane]
    }), axis=1)

    diff = pd.merge(t1, t2, on='S')  # Only get data for AT_CENTER
    if std:
        ax.plot(diff['S'], (diff['std_bl1']-diff['std_bl2'])/diff['std_bl1'], '^-', color=palette[plane],
                markeredgecolor=palette[plane], markersize=2, linewidth=1)
    if mean:
        ax.plot(diff['S'], 100*(diff['mean_bl1']-diff['mean_bl2']), '^-', color=palette[plane],
                markeredgecolor=palette[plane], markersize=2, linewidth=1)

    return diff