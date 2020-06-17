from .common import filled_plot
from .common import palette as common_palette
from .beam import halo_from_beam_o_df as compute_halo
import pandas as pd
import numpy as np


def tracking(ax, bl,beam_o_df, mean=True, std=True, halo=True, **kwargs):
    """Plot the beam envelopes from tracking data."""
    if kwargs.get("plane") is None:
        raise Exception("Plane (plane='X' or plane='Y') must be specified.")

    plane = kwargs.get("plane")
    palette = kwargs.get("palette", common_palette)
    if plane is None:
        raise Exception("The 'plane' keyword argument must be set to 'X' or 'Y'.")
    halo_99 = kwargs.get("halo99", halo)
    std_bpm = kwargs.get("std_bpm", False)

    dico_plane = {'X': 0, 'PX': 1, 'Y': 2, 'PY': 3}

    bl = bl.reset_index()

    t = bl.apply(lambda r: pd.Series({
        'NAME': r['NAME'],
        'S': r[kwargs.get("reference_plane", 'AT_EXIT')],
        '1%': 1000 * compute_halo(beam_o_df, r['NAME'], 0.01, plane) if halo_99 else 0.0,
        '5%': 1000 * compute_halo(beam_o_df, r['NAME'], 0.05, plane) if halo else 0.0,
        '95%': 1000 * compute_halo(beam_o_df, r['NAME'], 0.95, plane) if halo else 0.0,
        '99%': 1000 * compute_halo(beam_o_df, r['NAME'], 0.99, plane) if halo_99 else 0.0,
        'mean': 1000 * beam_o_df['BEAM_OUT'][r['NAME']][:, dico_plane[plane]].mean() if mean else 0.0,
        'std': 1000 * beam_o_df['BEAM_OUT'][r['NAME']][:, dico_plane[plane]].std() if std else 0.0,
        #'std_bpm': 1000 * r['BEAM'].std_bpm[plane][0] * int(pd.notnull(r['BPM'])) if 'BPM' in bl.line.columns and std_bpm else 0.0,
        #'std_bpm_err': np.max([1.0, 1000 * r['BEAM'].std_bpm[plane][1] * int( pd.notnull(r['BPM'])) if 'BPM' in bl.line.columns and std_bpm else 0.0]),
    }), axis=1)

    t.set_index('NAME')

    if t['S'].count == 0:
        return

    if halo:
        filled_plot(ax, t['S'], t['5%'], t['95%'], palette[plane], True, alpha=0.3)
        filled_plot(ax, t['S'], t['mean'] - t['std'], t['mean'] + t['std'], palette[plane], True, alpha=0.3)
        if halo_99:
            filled_plot(ax, t['S'], t['1%'], t['99%'], palette[plane], True, alpha=0.3)

    if std:
        ax.plot(t['S'], t['mean'] + t['std'],
                '^-',
                color=palette[plane],
                markeredgecolor=palette[plane],
                markersize=2,
                linewidth=1
                )
        ax.plot(t['S'], t['mean'] - t['std'],
                'v-',
                color=palette[plane],
                markeredgecolor=palette[plane],
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
                    color=palette['green'])
        ax.errorbar(t['S'] - 0.05, -t['std_bpm'], xerr=0.1, yerr=t['std_bpm_err'],
                    fmt='none',
                    elinewidth=2.0,
                    linewidth=0.0,
                    color=palette['green'])

    if mean:
        ax.plot(t['S'], t['mean'],
                '*-',
                color=palette[plane],
                markeredgecolor=palette[plane],
                markersize=2,
                linewidth=1,
                label=kwargs.get("label")
                )
