import numpy as np
from .constants import *


def aperture_check(b, e):
    # Circular aperture
    if e[INDEX_APERTYPE_CODE] == APERTYPE_CODE_CIRCLE:
        s = (b[:, 0]**2 + b[:, 2]**2) < e[INDEX_APERTURE]**2
    # Rectangular aperture
    elif e[INDEX_APERTYPE_CODE] == APERTYPE_CODE_RECTANGLE:
        s = (b[:, 0]**2 < e[INDEX_APERTURE]**2) & (b[:, 2]**2 < e[INDEX_APERTURE_2]**2)
    # Unknown aperture type
    else:
        return b
    # np.compress used for performance
    return np.compress(s, b, axis=0)