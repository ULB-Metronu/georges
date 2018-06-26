import numpy as np
from .constants import *


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

    t = np.zeros(4, 4, 4)
    t[111] = -(h/2)*np.tan(e1)**2
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


tensors = {
    CLASS_CODES['DRIFT']: None,
    CLASS_CODES['COLLIMATOR']: None,
    CLASS_CODES['SBEND']: sbend,
    CLASS_CODES['RBEND']: rbend,
    CLASS_CODES['QUADRUPOLE']: None,
    CLASS_CODES['ROTATION']: None,
}
