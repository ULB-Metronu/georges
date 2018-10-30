import numpy as np
from .constants import *


def sextupole(e):
    """
    Second order tensor of a pure sextupole
    :param e:
    :return:
    """
    k2 = e[INDEX_K2]
    k22 = k2**2
    length = e[INDEX_LENGTH]
    length2 = length**2
    length3 = length**3
    length4 = length**4
    t = np.zeros(5, 5, 5)
    t[0, 0, 0] = - 0.5 * k22 * length2
    t[0, 0, 1] = - (1.0/3.0) * k22 * length3
    t[0, 1, 1] = - (1.0/12.0) * k22 * length4
    t[0, 2, 2] = 0.5 * k22 * length2
    t[0, 2, 3] = (1.0/3.0) * k22 * length3
    t[0, 3, 3] = (1.0/12.0) * k22 * length4
    t[1, 0, 0] = - k22 * length
    t[1, 0, 1] = - k22 * length2
    t[1, 1, 1] = - (1.0/3.0) * k22 * length3
    t[1, 2, 2] = k22 * length
    t[1, 2, 3] = k22 * length2
    t[1, 3, 3] = (1.0/3.0) * k22 * length3
    t[2, 0, 2] = k22 * length2
    t[2, 0, 3] = (1.0/3.0) * k22 * length3
    t[2, 1, 2] = 0

def sbend(e):
    """
    Transfer matrix of a drift.
    :param e: element definition
    :return: a numpy array representing the 5D transfer matrix
    """
    theta = e[INDEX_ANGLE]
    length = e[INDEX_LENGTH]
    e1 = e[INDEX_E1]
    e2 = e[INDEX_E2]
    h1 = e[INDEX_H1]
    h2 = e[INDEX_H2]
    h = e[INDEX_ANGLE] / e[INDEX_LENGTH]

    t = np.zeros(5, 5, 5)
    t[0, 0, 0] = -(h/2)*np.tan(e1)**2
    t[133] = (h/2)*(1/np.cos(e1))**2
    t[211] = 0
    t[212] = 0
    t[233] = 0
    t[234] = 0
    t[313] = 0
    t[413] = 0
    t[414] = 0
    t[423] = 0
    return t


def rbend(e):
    """
    Tensor for a rbend element.
    :param e: element definition
    :return: a numpy array representing the element's tensor
    """
    t = np.zeros(4, 4, 4)
    return t


tensors = {
    CLASS_CODES['DRIFT']: None,
    CLASS_CODES['COLLIMATOR']: None,
    CLASS_CODES['SBEND']: sbend,
    CLASS_CODES['RBEND']: rbend,
    CLASS_CODES['QUADRUPOLE']: None,
    CLASS_CODES['ROTATION']: None,
}

