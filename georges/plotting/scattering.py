import matplotlib
import matplotlib.patches


def draw_slab(ax, e):
    materials_colors = {
        'graphite': 'g',
        'beryllium': 'r',
        'water': 'b',
    }
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], -1),  # (x,y)
            e['LENGTH'],  # width
            2,  # height
            hatch='', facecolor=materials_colors[e['MATERIAL']]
        )
    )


def draw_measuring_plane(ax, e):
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY']-0.005, -1),  # (x,y)
            0.01,  # width
            2,  # height
            hatch='', facecolor='k'
        )
    )


def scattering(ax, bl, **kwargs):
    bl.line.query("TYPE == 'slab'").apply(lambda e: draw_slab(ax, e), axis=1)
    bl.line.query("TYPE == 'mp'").apply(lambda e: draw_measuring_plane(ax, e), axis=1)
