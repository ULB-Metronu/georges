from .common import palette, filled_plot
import numpy as np
import pandas as pd


def tracking(ax, bl, **kwargs):
    """Plot the beam envelopes from tracking data."""
    if kwargs.get("plane") is None:
        raise Exception("Plane (plane='X' or plane='Y') must be specified.")

    plane = kwargs.get("plane")
    halo = kwargs.get("halo", True)
    std = kwargs.get("std", False)
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
    }), axis=1)

    if halo:
        filled_plot(ax, t['S'], t['1%'], t['99%'],palette[plane], True, alpha=0.3)
        filled_plot(ax, t['S'], t['5%'], t['95%'],palette[plane], True, alpha=0.3)
        filled_plot(ax,t['S'], -t['std'], t['std'],palette[plane], True, alpha=0.3)

    if std:
      ax.plot(t['S'], t['std'], '^-', color=palette[plane], markeredgecolor=palette[plane], linewidth=1.0)
      ax.plot(t['S'], -t['std'], 'v-', color=palette[plane], markeredgecolor=palette[plane], linewidth=1.0)

     if mean:
          ax.plot(t['S'], t['mean'], '*-', color=palette[plane], markeredgecolor=palette[plane], linewidth=1.0)
