from .common import palette, filled_plot
import numpy as np
import pandas as pd


def g4blprofile(ax, g4beamlineData, **kwargs):

    if kwargs.get("plane") is None:
        raise Exception("Plane (plane='X' or plane='Y') must be specified.")

    plane = kwargs.get("plane")
    std = kwargs.get("std", True)

    if plane == 'X':
        if std:
            # ax.plot(g4beamlineData['S'], g4beamlineData['sigmax'], color=palette[plane], markeredgecolor=palette[plane], linewidth=1.0)
            # ax.plot(g4beamlineData['S'], -g4beamlineData['sigmax'], color=palette[plane], markeredgecolor=palette[plane], linewidth=1.0)
            ax.plot(g4beamlineData['S'], g4beamlineData['sigmax'], '--k', linewidth=1.0)
            ax.plot(g4beamlineData['S'], -g4beamlineData['sigmax'], '--k', linewidth=1.0)
        # Mean
        ax.plot(g4beamlineData['S'], g4beamlineData['xmean'], '*-k', linewidth=1.0,label='G4Beamline')

    if plane == 'Y':
        if std:
            ax.plot(g4beamlineData['S'], g4beamlineData['sigmay'], '--k',  linewidth=1.0)
            ax.plot(g4beamlineData['S'], -g4beamlineData['sigmay'], '--k', linewidth=1.0)

        # Mean
        ax.plot(g4beamlineData['S'], g4beamlineData['ymean'], '*-k', linewidth=1.0,label='G4Beamline')


def plotg4enveloppe(ax, DataPlot):
    """ plot the enveloppe wich is defined by E(z)=eps*beta(z)"""
    # DataPlot[0]=mean
    # DataPlot[1]=eps
    # DataPlot[2]=beta

    DataPlot['Product']=np.sqrt(DataPlot['Emittance']*DataPlot['Beta'])
    enveloppe_Min=DataPlot['meanPos']-DataPlot['Product']
    enveloppe_Max=DataPlot['meanPos']+DataPlot['Product']

    ax.fill_between(DataPlot.index, DataPlot['meanPos']-enveloppe_Min, DataPlot['meanPos']+enveloppe_Max,color='blue', lw=1, alpha=0.5)
