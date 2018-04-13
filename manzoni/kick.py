import numpy as np
from .constants import *

try:
    import numpy.random_intel as nprandom
except ModuleNotFoundError:
    import numpy.random as nprandom


def sextupole(e, b, **kwargs):
    b[:, 1] = e[INDEX_K2] * (b[:, X]**2 - b[:, Y]**2)
    b[:, 3] = e[INDEX_K2] * (b[:, X] * b[:, Y])
    return b


def octupole(e, b, **kwargs):
    b[:, 1] = e[INDEX_K3] * (b[:, X]**3 - 3 * b[:, X] * b[:, Y]**2)
    b[:, 3] = e[INDEX_K3] * (3 * b[:, X]**2 * b[:, Y] - b[:, Y]**3)
    return b



def decapole(e, b, **kwargs):
    b[:, 1] = e[INDEX_K3] * (b[:, X]**4 - 6 * b[:, X]**2 * b[:, Y]**2 + b[:, Y]**4)
    b[:, 3] = e[INDEX_K3] * (b[:, X]**3 * b[:, Y] - b[:, X] * b[:, Y]**3)
    return b


def multipole(e, b, **kwargs):
    return np.array(
        [
            0,
            0,
            0,
            0,
            0
        ]
    )


def hkicker(e, b, **kwargs):
    return b + np.array(
        [
            0,
            e[INDEX_K1],
            0,
            0,
            0
        ]
    )


def vkicker(e, b, **kwargs):
    return b + np.array(
        [
            0,
            0,
            0,
            e[INDEX_K1],
            0
        ]
    )


def degrader(e, b, **kwargs):
    # Remove particles
    idx = np.random.randint(b.shape[0], size=int((1 - e[INDEX_FE_LOSS]) * b.shape[0]))
    b = b[idx, :]
    b += nprandom.multivariate_normal(
        [0.0, 0.0, 0.0, 0.0, 0.0],
        np.array(
            [
                [e[INDEX_FE_A0], e[INDEX_FE_A1], 0, 0, 0],
                [e[INDEX_FE_A1], e[INDEX_FE_A2], 0, 0, 0],
                [0, 0, e[INDEX_FE_A0], e[INDEX_FE_A1], 0],
                [0, 0, e[INDEX_FE_A1], e[INDEX_FE_A2], 0],
                [0, 0, 0, 0, e[INDEX_FE_DPP]]
            ]
        ),
        int(b.shape[0]))
    return b


kick = {
    CLASS_CODES['SEXTUPOLE']: sextupole,
    CLASS_CODES['OCTUPOLE']: octupole,
    CLASS_CODES['DECAPOLE']: decapole,
    CLASS_CODES['MULTIPOLE']: multipole,
    CLASS_CODES['HKICKER']: hkicker,
    CLASS_CODES['VKICKER']: vkicker,
    CLASS_CODES['DEGRADER']: degrader,
}
