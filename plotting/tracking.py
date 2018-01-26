from .common import palette, filled_plot
import pandas as pd
import numpy as np


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
        'std_bpm': 1000 * r['BEAM'].std_bpm[plane][0] if 'BPM' in bl.line.columns else -10000.0,
        'std_bpm_err': np.max([1.0, 1000 * r['BEAM'].std_bpm[plane][1] if 'BPM' in bl.line.columns else 0.0]),
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
        ax.errorbar(t['S']- 0.05, t['std_bpm'], xerr=0.1, yerr=t['std_bpm_err'],
                    fmt='none',
                    elinewidth=2.0,
                    linewidth=0.0,
                    color=palette['green'])
        ax.errorbar(t['S']- 0.05, -t['std_bpm'], xerr=0.1, yerr=t['std_bpm_err'],
                    fmt='none',
                    elinewidth=2.0,
                    linewidth=0.0,
                    color=palette['green'])

    if mean:
        ax.plot(t['S'], t['mean'], '*-', color=palette[plane],
                markeredgecolor=palette[plane], markersize=2, linewidth=1, label=kwargs.get("label"))


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
        ax.plot(g4beamlineData['S'], g4beamlineData['xmean'], '*-k', linewidth=1.0, label='G4Beamline')

    if plane == 'Y':
        if std:
            ax.plot(g4beamlineData['S'], g4beamlineData['sigmay'], '--k', linewidth=1.0)
            ax.plot(g4beamlineData['S'], -g4beamlineData['sigmay'], '--k', linewidth=1.0)

        # Mean
        ax.plot(g4beamlineData['S'], g4beamlineData['ymean'], '*-k', linewidth=1.0, label='G4Beamline')


def plotg4enveloppe(ax, DataPlot):
    """ plot the enveloppe wich is defined by E(z)=eps*beta(z)"""
    # DataPlot[0]=mean
    # DataPlot[1]=eps
    # DataPlot[2]=beta

    DataPlot['Product'] = np.sqrt(DataPlot['Emittance'] * DataPlot['Beta'])
    enveloppe_Min = DataPlot['meanPos'] - DataPlot['Product']
    enveloppe_Max = DataPlot['meanPos'] + DataPlot['Product']

    ax.fill_between(DataPlot.index, DataPlot['meanPos'] - enveloppe_Min, DataPlot['meanPos'] + enveloppe_Max,
                    color='blue', lw=1, alpha=0.5)
