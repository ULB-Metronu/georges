from plotting.common import palette, filled_plot
from beam import Beam
from plotting.common import palette, filled_plot
import numpy as np


def track(ax, bl, plane):
    bl = bl.line
    filled_plot(ax, bl.index, 1000 * bl['BEAM'].apply(lambda x: x.mean())['X'] + 1000 * bl['BEAM'].apply(lambda x: x.std())['X'])


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
			 
def plotg4enveloppe(ax,DataPlot):
    """ plot the enveloppe wich is defined by E(z)=eps*beta(z)"""
    #DataPlot[0]=mean
    #DataPlot[1]=eps
    #DataPlot[2]=beta

    DataPlot['Product']=np.sqrt(DataPlot['Emittance']*DataPlot['Beta'])
    enveloppe_Min=DataPlot['meanPos']-DataPlot['Product']
    enveloppe_Max=DataPlot['meanPos']+DataPlot['Product']
    
    ax.fill_between(DataPlot.index, DataPlot['meanPos']-enveloppe_Min, DataPlot['meanPos']+enveloppe_Max,color='blue', lw=1, alpha=0.5)

