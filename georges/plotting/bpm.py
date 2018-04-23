from .common import palette
import pandas as pd
import matplotlib.patches as patches


def draw_bpm_size(ax, s, x):
    ax.add_patch(
        patches.Rectangle(
            (s - 0.05, -x),
            0.1,
            2 * x,
        )
    )


def bpm(ax, bl, **kwargs):
    """TODO."""
    if kwargs.get('plane') is None:
        raise Exception("'plane' argument must be provided.")
    bl.line[bl.line[f"BPM_STD_{kwargs.get('plane')}"].notnull()].apply(
        lambda x: draw_bpm_size(ax, x['AT_CENTER'], x[f"BPM_STD_{kwargs.get('plane')}"]),
        axis=1
    )
