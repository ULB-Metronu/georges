import matplotlib
import matplotlib.patches
import numpy as np
from .common import palette

STYLE_BEND_HATCH = '/'
STYLE_QUAD_HATCH = ''


def xy_from_string(a, i, c):
    def convert(x):
        try:
            converted = float(x)
        except ValueError:
            try:
                converted = float(c.get(x))
            except TypeError:
                converted = 1.0
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
            e['L'].magnitude,  # width
            1000 * e['CHAMBER_UP'],  # height
            hatch='', facecolor=palette['base01']
        )
    )
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], 1000 * (-e['APERTURE_DOWN'])),  # (x,y)
            e['L'].magnitude,  # width
            -1000 * e['CHAMBER_UP'],  # height
            hatch='', facecolor=palette['base01']
        )
    )


def draw_quad(ax, e):
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], 1000 * e['APERTURE_UP'] + e['CHAMBER_UP']),  # (x,y)
            e['L'].magnitude,  # width
            100,  # height
            hatch=STYLE_QUAD_HATCH, facecolor=palette['quad']
        )
    )

    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], -1000 * e['APERTURE_DOWN'] - e['CHAMBER_UP']),  # (x,y)
            e['L'].magnitude,  # width
            -100,  # height
            hatch=STYLE_QUAD_HATCH, facecolor=palette['quad']
        )
    )
    draw_chamber(ax, e)


def draw_coll(ax, e, plane):
    #if 'PIPE' not in e:
    #    return
    #if not np.isnan(e['PIPE']):
    #    return
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], 1000 * e['APERTURE_UP']),  # (x,y)
            e['L'].magnitude,  # width
            100,  # height
            facecolor=palette['coll']
        )
    )

    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], -1000 * e['APERTURE_DOWN']),  # (x,y)
            e['L'].magnitude,  # width
            -100,  # height
            facecolor=palette['coll']
        )
    )


def draw_bend(ax, e):
    tmp = 1000 * e['APERTURE_UP'] + e['CHAMBER_UP']
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], tmp if tmp < 55 else 55),  # (x,y)
            e['L'].magnitude,  # width
            100,  # height
            hatch=STYLE_BEND_HATCH, facecolor=palette['bend']
        )
    )
    tmp = -1000 * e['APERTURE_DOWN'] - e['CHAMBER_UP']
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (e['AT_ENTRY'], tmp if abs(tmp) < 55 else -55),  # (x,y)
            e['L'].magnitude,  # width
            -100,  # height
            hatch=STYLE_BEND_HATCH, facecolor=palette['bend']
        )
    )
    draw_chamber(ax, e)


def fill_aperture(element, context):
    if element.name+'_APERTURE' in context and element['TYPE'] == 'SLITS':
        element['APERTURE'] = context[element.name + '_APERTURE']
    if element['TYPE'] == 'COLLIMATOR' and \
            element['PLUG'] == 'APERTURE' and \
            element['APERTURE'] is not None and \
            np.isnan(element['APERTURE']):
        element['APERTURE'] = element['CIRCUIT']
    return element


def aperture(ax, bl, **kwargs):

    if 'APERTURE' not in bl:
        return

    planes = kwargs.get('plane', 'both')

    if planes == 'X':
        index = 0
    elif planes == 'Y':
        index = 1

    bl['APERTURE_UP'] = bl['APERTURE'].apply(
        lambda a: a[index].magnitude
    )
    bl['APERTURE_DOWN'] = bl['APERTURE'].apply(
        lambda a: a[index].magnitude
    )

    if 'CHAMBER' not in bl:
        bl['CHAMBER'] = 0

    bl['CHAMBER_UP'] = bl['CHAMBER'].apply(
        lambda a: a
    )
    bl['CHAMBER_DOWN'] = bl['CHAMBER'].apply(
        lambda a: a
    )

    bl.query("CLASS == 'Quadrupole'").apply(lambda e: draw_quad(ax, e), axis=1)
    bl.query("CLASS == 'SBend'").apply(lambda e: draw_bend(ax, e), axis=1)
    bl.query("CLASS == 'RBend'").apply(lambda e: draw_bend(ax, e), axis=1)
    bl.query("CLASS == 'RectangularCollimator'").apply(lambda e: draw_coll(ax, e, planes), axis=1)
    bl.query("CLASS == 'CircularCollimator'").apply(lambda e: draw_coll(ax, e, planes), axis=1)