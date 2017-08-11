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
        's11': s11,
        's12': s12,
        's22': s22,
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
