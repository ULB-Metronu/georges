from .common import palette, filled_plot
import pandas as pd


def tracking(ax, bl, **kwargs):
    """Plot the beam envelopes from tracking data."""
    if kwargs.get("plane") is None:
        raise Exception("Plane (plane='X' or plane='Y') must be specified.")

    plane = kwargs.get("plane")
    if plane is None:
        raise Exception("The 'plane' keyword argument must be set to 'X' or 'Y'.")
    halo = kwargs.get("halo", True)
    halo_99 = kwargs.get("halo99", halo)
    std = kwargs.get("std", True)
    std_bpm = kwargs.get("std_bpm", False)
    mean = kwargs.get("mean", True)

    t = bl.line.query("BEAM == BEAM").apply(lambda r: pd.Series({
        'S': r['AT_CENTER'],
        '1%': 1000 * r['BEAM'].halo['1%'][plane],
        '5%': 1000 * r['BEAM'].halo['5%'][plane],
        '20%': 1000 * r['BEAM'].halo['20%'][plane],
        '80%': 1000 * r['BEAM'].halo['80%'][plane],
        '95%': 1000 * r['BEAM'].halo['95%'][plane],
        '99%': 1000 * r['BEAM'].halo['99%'][plane],
        'mean': 1000 * r['BEAM'].mean[plane],
        'std': 1000 * r['BEAM'].std[plane],
        'std_bpm': 1000 * r['BEAM'].std_bpm[plane] if r['CLASS'] == 'MARKER' else 0.0,
    }), axis=1)

    if t['S'].count == 0:
        return

    if halo:
        if halo_99:
            filled_plot(ax, t['S'], t['1%'], t['99%'], palette[plane], True, alpha=0.3)
        filled_plot(ax, t['S'], t['5%'], t['95%'], palette[plane], True, alpha=0.3)
        filled_plot(ax, t['S'], -t['std'], t['std'], palette[plane], True, alpha=0.3)

    if std:
        ax.plot(t['S'], t['std'], '^-', color=palette[plane],
                markeredgecolor=palette[plane], markersize=2, linewidth=1)
        ax.plot(t['S'], -t['std'], 'v-', color=palette[plane],
                markeredgecolor=palette[plane], markersize=2, linewidth=1)

    if std_bpm:
        ax.plot(t['S'], t['std_bpm'], '^', color=palette['green'],
                markeredgecolor=palette['green'], markersize=2, linewidth=1)
        ax.plot(t['S'], -t['std_bpm'], 'v', color=palette['green'],
                markeredgecolor=palette['green'], markersize=2, linewidth=1)

    if mean:
        ax.plot(t['S'], t['mean'], '*-', color=palette[plane],
                markeredgecolor=palette[plane], markersize=2, linewidth=1, label=kwargs.get("label"))
