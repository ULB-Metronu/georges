from plotting.common import palette, filled_plot
from georges.beam import Beam


def track(ax, bl, plane):
    bl = bl.line
    print(bl['BEAM'])
    filled_plot(ax, bl.index, 1000 * bl['BEAM'].apply(lambda x: x.std()))
    pass


def tracking(ax, bl, context, **kwargs):
    """Plot the beam envelopes from tracking data."""
    if kwargs.get("plane") is None:
        raise Exception("Plane (plane='X' or plane='Y') must be specific.")


    filled_plot(ax, envelope.index, 1000 * trajectory[plane], 1000 * trajectory[plane] + 1000 * envelope[plane],
                       palette[plane], True, alpha=0.4)
    filled_plot(ax, envelope.index, 1000 * trajectory[plane], 1000 * trajectory[plane] - 1000 * envelope[plane],
                       palette[plane], True, alpha=0.4)

    filled_plot(ax, envelope.index, 1000 * trajectory[plane], 1000 * halo_sup[plane],
                       palette[plane], True, alpha=0.2)
    filled_plot(ax, envelope.index, 1000 * trajectory[plane], 1000 * halo_inf[plane],
                       palette[plane], True, alpha=0.2)
    filled_plot(ax, envelope.index, 1000 * trajectory[plane], 1000 * halo_sup_bis[plane],
                       palette[plane], True, alpha=0.1)
    filled_plot(ax, envelope.index, 1000 * trajectory[plane], 1000 * halo_inf_bis   [plane],
                       palette[plane], True, alpha=0.1)

    ax.plot(envelope2.index, 1000 * trajectory[plane] + 1000 * envelope[plane], '*-',
             color=palette[plane],
             markeredgecolor=palette[plane],
             linewidth=1.0, alpha=0.4)
    ax.plot(envelope2.index, 1000 * trajectory[plane] - 1000 * envelope[plane], '*-',
             color=palette[plane],
             markeredgecolor=palette[plane],
             linewidth=1.0, alpha=0.4)

    ax.plot(trajectory.index, 1000 * trajectory[plane], '*-',
             color=palette[plane],
             markeredgecolor=palette[plane],
             linewidth=1.0)