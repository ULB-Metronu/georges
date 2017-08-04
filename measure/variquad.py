import numpy as np
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


def variquad(**kargs):
    """Implementation of the quadrupole scan 'variquad' method"""
    with_plot = kwargs.get("plot", None)
    start = kwargs.get("start")
    end = kwargs.get("end")
    plane = kwargs.get("plane", 'X')
    quad_length = kwargs.get("quad_length", 1.0)
    bl = kwargs.get("line")
    context = kwargs.get("context")
    data = kwargs.get("data")
    if data is None:
        raise VariquadException("data must be provided as a numpy.array")

    x = variquad[:,0]
    y = variquad[:,1]**2
    popt, pcov = variquad_fit(x, y)
    if with_plot is not None:
        with_plot.plot(x, y**2, '*', label='data')
        with_plot.plot(x, quadratic(x, *popt), '-', label='fit')
        with_plot.legend()

    bl_map = georges.madx.sectormap(line=bl, context=context, start=start, places=[start, end])
    if plane == 'X':
        r11 = bl_map.loc[end]['R33']
        r12 = bl_map.loc[end]['R34']
        r21 = bl_map.loc[end]['R43']
        r22 = bl_map.loc[end]['R44']
    elif plane == 'Y':
        r11 = bl_map.loc[end]['R33']
        r12 = bl_map.loc[end]['R34']
        r21 = bl_map.loc[end]['R43']
        r22 = bl_map.loc[end]['R44']
    else:
        raise VariquadException("plane must be 'X' or 'Y'")

    s11 = np.sqrt(popt[0] / (quad_length * (r12**2)))
    s12 = 1
    s22 = 1
    eps = 1

    return {
        's11': s11,
        's12': s12,
        's22': s22,
        'eps': eps,
        'fit': [popt, pcov],
        'map': bl_map
    }
