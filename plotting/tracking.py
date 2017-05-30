from georges.plotting.common import palette, filled_plot


def track(ax, bl, plane):
    bl = bl.line
    print(bl['BEAM'])
    filled_plot(ax, bl.index, 1000 * bl['BEAM'].apply(lambda x: x.std()))
    pass


def tracking(ax, tracking, plane):
    envelope = tracking['envelope']
    envelope2 = tracking['envelope2']
    trajectory = tracking['trajectory']
    halo_sup = tracking['halo_sup']
    halo_inf = tracking['halo_inf']
    halo_sup_bis = tracking['halo_sup_bis']
    halo_inf_bis = tracking['halo_inf_bis']
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
			 
def g4tracking(ax, tracking, plane):

    envelope = tracking['envelope']
    envelope2 = tracking['envelope2']
    trajectory = tracking['trajectory']
    halo_sup = tracking['halo_sup']
    halo_inf = tracking['halo_inf']
    halo_sup_bis = tracking['halo_sup_bis']
    halo_inf_bis = tracking['halo_inf_bis']
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