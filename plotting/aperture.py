import matplotlib
import matplotlib.patches
import numpy as np
from .common import palette


def xy_from_string(a, i, c):
    def convert(x):
        try:
            converted = float(x)
        except ValueError:
            converted = float(c.get(x))
        return converted
    if len(str(a).strip('[]').strip('{}').split(',')) >= int(i) + 1:
        return convert(str(a).strip('[]').strip('{}').split(',')[int(i)])
    elif len(str(a).strip('[]').strip('{}').split(',')) > 0:
        return convert(str(a).strip('[]').strip('{}').split(',')[0])
    else:
        return np.inf


def draw_chamber(ax, e):
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], 1000 * (e['APERTURE_UP'])),  # (x,y)
            e['ORBIT_LENGTH'],  # width
            1000 * e['CHAMBER_UP'],  # height
            hatch='', facecolor=palette['base01']
        )
    )
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], 1000 * (-e['APERTURE_DOWN'])),  # (x,y)
            e['ORBIT_LENGTH'],  # width
            -1000 * e['CHAMBER_UP'],  # height
            hatch='', facecolor=palette['base01']
        )
    )


def draw_quad(ax, e):
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], 1000 * e['APERTURE_UP'] + e['CHAMBER_UP']),  # (x,y)
            e['ORBIT_LENGTH'],  # width
            100,  # height
            hatch='.', facecolor=palette['quad']
        )
    )

    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], -1000 * e['APERTURE_DOWN'] - e['CHAMBER_UP']),  # (x,y)
            e['ORBIT_LENGTH'],  # width
            -100,  # height
            hatch='.', facecolor=palette['quad']
        )
    )
    draw_chamber(ax, e)


def draw_coll(ax, e, plane):
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], 1000 * e['APERTURE_UP']),  # (x,y)
            e['LENGTH'],  # width
            100,  # height
            facecolor=palette['coll']
        )
    )

    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], -1000 * e['APERTURE_DOWN']),  # (x,y)
            e['LENGTH'],  # width
            -100,  # height
            facecolor=palette['coll']
        )
    )


def draw_bend(ax, e):
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], 1000 * e['APERTURE_UP'] + e['CHAMBER_UP']),  # (x,y)
            e['ORBIT_LENGTH'],  # width
            100,  # height
            hatch='/', facecolor=palette['bend']
        )
    )
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], -1000 * e['APERTURE_DOWN'] - e['CHAMBER_UP']),  # (x,y)
            e['ORBIT_LENGTH'],  # width
            -100,  # height
            hatch='/', facecolor=palette['bend']
        )
    )
    draw_chamber(ax, e)


def fill_aperture(element, context):
    if element.name+'_APERTURE' in context and element['TYPE'] == 'SLITS':
        element['APERTURE'] = context[element.name+'_APERTURE']
    if element['TYPE'] == 'COLLIMATOR' and element['PLUG'] == 'APERTURE':
        element['APERTURE'] = element['CIRCUIT']
    return element


def aperture(ax, bl, **kwargs):
    if kwargs.get('context'):
        bl = bl.line.apply(lambda e: fill_aperture(e, kwargs.get('context')), axis=1)
    else:
        bl = bl.line

    if 'APERTURE' not in bl:
        return

    planes = kwargs.get('plane', 'both')

    bl['APERTURE_UP'] = bl['APERTURE'].apply(
        lambda a: xy_from_string(a, planes == 'both' or planes == 'Y', kwargs.get('context'))
    )
    bl['APERTURE_DOWN'] = bl['APERTURE'].apply(
        lambda a: xy_from_string(a, not (planes == 'both' or planes == 'X'), kwargs.get('context'))
    )

    if 'CHAMBER' not in bl:
        bl['CHAMBER'] = 0

    bl['CHAMBER_UP'] = bl['CHAMBER'].apply(
        lambda a: xy_from_string(a, planes == 'both' or planes == 'Y', kwargs.get('context'))
    )
    bl['CHAMBER_DOWN'] = bl['CHAMBER'].apply(
        lambda a: xy_from_string(a, not (planes == 'both' or planes == 'X'), kwargs.get('context'))
    )

    bl.query("CLASS == 'QUADRUPOLE'").apply(lambda e: draw_quad(ax, e), axis=1)
    bl.query("CLASS == 'SBEND'").apply(lambda e: draw_bend(ax, e), axis=1)
    bl.query("CLASS == 'RBEND'").apply(lambda e: draw_bend(ax, e), axis=1)
    bl.query("CLASS == 'COLLIMATOR'").apply(lambda e: draw_coll(ax, e, planes), axis=1)
