from .common import palette, filled_plot
import numpy as np
import pandas as pd


def tracking(ax, bl, context, **kwargs):
    """Plot the beam envelopes from tracking data."""
    if kwargs.get("plane") is None:
        raise Exception("Plane (plane='X' or plane='Y') must be specified.")

    plane = kwargs.get("plane")
    bl.line=(bl.line[bl.line['TYPE'] == 'MARKER']) # To remove the NAN in beam
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

    filled_plot(ax, t['S'], t['1%'], t['99%'],palette[plane], True, alpha=0.3)
    filled_plot(ax, t['S'], t['5%'], t['95%'],palette[plane], True, alpha=0.3)
    filled_plot(ax, t['S'], t['20%'], t['80%'],palette[plane], True, alpha=0.3)

    ax.plot(t['S'], t['mean'], '*-', color=palette[plane], markeredgecolor=palette[plane], linewidth=1.0)

    # if(plane =='X'):
    #
    #     ax.plot(t['S'], -t['std'], '*-', color=palette[plane], markeredgecolor=palette[plane], linewidth=1.0)
    # else:
    #     ax.plot(t['S'], t['std'], '*-', color=palette[plane], markeredgecolor=palette[plane], linewidth=1.0)




    #ax.plot(trajectory.index, 1000 * trajectory[plane], '*-',
     #        color=palette[plane],
      #       markeredgecolor=palette[plane],
       #      linewidth=1.0)


def plotg4enveloppe(ax, DataPlot):
    """ plot the enveloppe wich is defined by E(z)=eps*beta(z)"""
    #DataPlot[0]=mean
    #DataPlot[1]=eps
    #DataPlot[2]=beta

    DataPlot['Product']=np.sqrt(DataPlot['Emittance']*DataPlot['Beta'])
    enveloppe_Min=DataPlot['meanPos']-DataPlot['Product']
    enveloppe_Max=DataPlot['meanPos']+DataPlot['Product']

    ax.fill_between(DataPlot.index, DataPlot['meanPos']-enveloppe_Min, DataPlot['meanPos']+enveloppe_Max,color='blue', lw=1, alpha=0.5)
