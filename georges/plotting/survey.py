import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd


def survey_iba(ax, bl, **kwargs):
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_xlim([0, np.max(bl.line['X']) / 1000])
    ax.set_ylim([0, np.max(bl.line['Y']) / 1000])
    for index, row in bl.line.iterrows():
        if pd.notnull(row['X']) and pd.notnull(row['Y']):
            if kwargs.get("labels", False):
                ax.annotate(index,
                            xy=(row['X'] / 1000, row['Y'] / 1000), xytext=(row['X'] / 1000 + 0.5, row['Y'] / 1000 + 0.5),
                            arrowprops=dict(arrowstyle="->", facecolor='black', shrinkA=50, shrinkB=5000, ),
                            size=3,
                            horizontalalignment='left',
                            verticalalignment='bottom',
                            clip_on=True
                            )
            if row['CLASS'] == 'QUADRUPOLE':
                ax.add_patch(
                    patches.Rectangle(
                        (row['X'] / 1000, row['Y'] / 1000),  # (x,y)
                        0.20,  # width
                        0.20,  # height
                        facecolor='#268bd2',
                        edgecolor='#268bd2'
                    )
                )
            elif row['CLASS'] == 'RBEND' or row['CLASS'] == 'SBEND':
                ax.add_patch(
                    patches.Rectangle(
                        (row['X'] / 1000, row['Y'] / 1000),  # (x,y)
                        0.250,  # width
                        0.250,  # height
                        facecolor='#dc322f',
                        edgecolor='#dc322f'
                    )
                )
            else:
                ax.add_patch(
                    patches.Rectangle(
                        (row['X'] / 1000, row['Y'] / 1000),  # (x,y)
                        0.250,  # width
                        0.250,  # height
                        facecolor='#657b83',
                        edgecolor='#657b83'
                    )
                )
    plt.axes().set_aspect('equal', 'datalim')


ELEMENT_COLORS = {
    'QUADRUPOLE': '#cb4b16',  # Solarized 'orange'
    'SEXTUPOLE': 'g',
    'OCTUPOLE': 'g',
    'MULTIPOLE': 'g'
}


def survey_madx(ax, bl):
    x_edges = [
        np.min(list(map(np.abs, [bl.line['Z'].max(), bl.line['Z'].min()]))),
        np.max(list(map(np.abs, [bl.line['Z'].max(), bl.line['Z'].min()])))
    ]
    y_edges = [
        np.min(list(map(np.abs, [bl.line['X'].max(), bl.line['X'].min()]))),
        np.max(list(map(np.abs, [bl.line['X'].max(), bl.line['X'].min()])))
    ]
    edges = [
        np.max([x_edges[0], y_edges[0]]),
        np.max([x_edges[1], y_edges[1]])
    ]
    padding = np.max(edges)*0.1
    ax.set_xlim([edges[0]-padding, edges[1]+padding])
    ax.set_ylim([edges[0]-padding, edges[1]+padding])
    tmp = -90
    for index, row in bl.line.iterrows():
        if row['KEYWORD'] == 'SBEND' or row['KEYWORD'] == 'MARKER':
            if row['KEYWORD'] == 'SBEND':
                r = row['L'] / row['ANGLE']
                if row.get('MAIN_BEND'):
                    c = 'r' if row['MAIN_BEND'] is True else 'b'
                else:
                    if row['ANGLE'] > 0:
                        c = 'r'
                    elif row['ANGLE'] < 0:
                        c = 'b'
                    else:
                        c = 'gray'
            if row['KEYWORD'] == 'MARKER':
                r = 0.1
                c = 'b'
            theta = -row['THETA']
            angle = -row['ANGLE'] / np.pi * 180.0
            centre = [(row['Z'] - np.sin(theta) * r), -(row['X'] - np.cos(theta) * r)]
            theta1 = min(tmp, tmp - angle)
            theta2 = max(tmp, tmp - angle)
            tmp = theta1 if angle > 0 else theta2
            w = patches.Wedge(centre, r + 0.2, theta1, theta2, width=0.4, alpha=1, facecolor=c, ec=c)
        if row['KEYWORD'] == 'DRIFT':
            w = patches.Rectangle(
                (
                    row['Z'] - row['L'] * np.cos(row['THETA']) - 0.2 * np.sin(row['THETA']),
                    -row['X'] + row['L'] * np.sin(row['THETA']) - 0.2 * np.cos(row['THETA'])
                ),
                row['L'],
                0.4,
                angle=np.degrees(-row['THETA']),
                alpha=0.2,
                facecolor='y',
                ec='y',
                hatch=''
            )
        if row['KEYWORD'] in ('MULTIPOLE', 'QUADRUPOLE', 'SEXTUPOLE', 'OCTUPOLE'):
            w = patches.Rectangle(
                (
                    row['Z'] - row['L'] * np.cos(row['THETA']) - 0.2 * np.sin(row['THETA']),
                    -row['X'] + row['L'] * np.sin(row['THETA']) - 0.2 * np.cos(row['THETA'])
                ),
                row['L'],
                0.4,
                angle=np.degrees(-row['THETA']),
                alpha=1.0,
                facecolor=ELEMENT_COLORS[row['KEYWORD']],
                ec=ELEMENT_COLORS[row['KEYWORD']],
                hatch=''
            )
        ax.add_patch(w)


def survey(ax, bl, style='madx', **kwargs):
    if style == 'iba_survey':
        survey_iba(ax, bl, **kwargs)
    elif style == 'madx':
        survey_madx(ax, bl)
    else:
        raise Exception("Style not supported.")
