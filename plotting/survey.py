import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd


def survey(ax, bl, **kwargs):
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
