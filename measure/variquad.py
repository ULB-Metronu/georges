from ..madx import sectormap
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


class VariquadException(Exception):
    """Exception raised for errors in the Variquad module."""

    def __init__(self, m):
        self.message = m


def quadratic(x, a, b, c):
    return a * x ** 2 + b * x + c


def variquad_fit(x, y):
    popt, pcov = curve_fit(quadratic, x, y, p0=[0, 0, 1])

    return [popt, np.sqrt(pcov.diagonal())]


def variquad(**kwargs):
    """Implementation of the quadrupole scan 'variquad' method"""
    start = kwargs.get("start")
    end = kwargs.get("end")
    plane = kwargs.get("plane", 'X')
    ql = kwargs.get("quad_length", 1.0)
    bl = kwargs.get("line")
    context = kwargs.get("context")
    debug = kwargs.get("debug", False)
    if kwargs.get("data") is None:
        raise VariquadException("data must be provided as a numpy.array")

    data = np.array([
        kwargs.get("data")[start].apply(lambda v: kwargs.get("i2k", lambda i: i)(v)).values,
        kwargs.get("data")["{}_fit_sigma_{}".format(end, plane)]/1000.0
    ]).transpose()
    x = data[:, 0]
    y = data[:, 1]**2
    popt, pcov = variquad_fit(x, y)

    bl_map = sectormap(line=bl, context=context, start=start, places=[start, end], debug=debug)
    if plane == 'X':
        r11 = bl_map.line.loc[end]['R11']
        r12 = bl_map.line.loc[end]['R12']
    elif plane == 'Y':
        r11 = bl_map.line.loc[end]['R33']
        r12 = bl_map.line.loc[end]['R34']
    else:
        raise VariquadException("plane must be 'X' or 'Y'")

    a = popt[0]
    b = popt[1]
    c = popt[2]
    s11 = a / (ql**2 * r12**2)
    s12 = (b - 2 * s11 * ql * r11 * r12)/(2*ql*r12**2)
    s21 = s12
    s22 = (c-s11 * r11**2 - 2 * s12 * r11 * r12)/(r12**2)
    emit = np.sqrt(np.linalg.det(np.array([[s11, s12], [s21, s22]])))
    beta = s11/emit
    alpha = -s12/emit
    gamma = s22/emit

    return {
        'S11': s11,
        'S12': s12,
        'S22': s22,
        'emit': emit,
        'beta': beta,
        'alpha': alpha,
        'gamma': gamma,
        'fit': [popt, pcov],
        'map': bl_map,
        'x': x,
        'y': y,
        'fit_function': (lambda u: quadratic(u, *popt))
    }


def variquad_plot(variquad_results, **kwargs):
    if kwargs.get("ax"):
        ax = kwargs.get("ax")
    else:
        ax = plt.figure().add_subplot(111)

    ax.plot(variquad_results['x'], variquad_results['y'], '*')
    ax.plot(variquad_results['x'], variquad_results['fit_function'](variquad_results['x']))


def backtrack(**kwargs):
    track_from = kwargs.get('track_from')
    track_to = kwargs.get('track_to')
    bl = kwargs.get("line")
    context = kwargs.get("context")
    plane = kwargs.get("plane")

    bl_map = sectormap(line=bl, context=context, start=track_to, places=[track_to, track_from])

    if plane == 'X':
        sigma_matrix = np.array([[context['S11'], context['S12']], [context['S12'], context['S22']]])
        tmp = bl_map.line.loc[track_from][['R11', 'R12', 'R21', 'R22']].values
    elif plane == 'Y':
        sigma_matrix = np.array([[context['S33'], context['S34']], [context['S34'], context['S44']]])
        tmp = bl_map.line.loc[track_from][['R33', 'R34', 'R43', 'R44']].values
    else:
        raise Exception("Invalid plane. 'plane' must be 'X' or 'Y'.")
    r_matrix = np.array([[tmp[0], tmp[1]], [tmp[2], tmp[3]]])
    inv_r_matrix = np.linalg.inv(r_matrix)
    res = np.matmul(inv_r_matrix, np.matmul(sigma_matrix, inv_r_matrix.transpose()))
    emit = np.sqrt(np.linalg.det(res))
    beta = res[0, 0] / emit
    alpha = -res[0, 1] / emit
    gamma = res[1, 1] / emit

    return {
        'emit': emit,
        'beta': beta,
        'alpha': alpha,
        'gamma': gamma,
        'S11': res[0, 0],
        'S12': res[0, 1],
        'S22': res[1, 1]
    }
