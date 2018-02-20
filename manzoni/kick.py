import numpy as np
from .constants import *


try:
    import numpy.random_intel as nprandom
except ModuleNotFoundError:
    import numpy.random as nprandom


def sextupole(e, b, **kwargs):
    return np.array(
        [
            0,
            0,
            0,
            0,
            0
        ]
    )


def octupole(e, b, **kwargs):
    return np.array(
        [
            0,
            0,
            0,
            0,
            0
        ]
    )


def decapole(e, b, **kwargs):
    return np.array(
        [
            0,
            0,
            0,
            0,
            0
        ]
    )


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
    return np.array(
        [
            0,
            e[INDEX_K1],
            0,
            0,
            0
        ]
    )


def vkicker(e, b, **kwargs):
    return np.array(
        [
            0,
            0,
            0,
            e[INDEX_K1],
            0
        ]
    )


def degrader(e, b, **kwargs):
    s11 = kwargs.get('deg', {}).get('A', [0, 0, 0])[0]
    s12 = kwargs.get('deg', {}).get('A', [0, 0, 0])[1]
    s22 = kwargs.get('deg', {}).get('A', [0, 0, 0])[2]
    dpp = kwargs.get('deg', {}).get('DPP', 0)

    # Remove particles
    idx = np.random.randint(b.shape[0], size=int((1 - kwargs.get('loss', 0)) * b.shape[0]))
    b = b[idx, :]

    return nprandom.multivariate_normal(
        [0.0, 0.0, 0.0, 0.0, 0.0], np.array(
            [
                [s11, s12, 0, 0, 0],
                [s12, s22, 0, 0, 0],
                [0, 0, s11, s12, 0],
                [0, 0, s12, s22, 0],
                [0, 0, 0, 0, dpp]
            ]),
        int(b.shape[0]))


kick = {
    CLASS_CODES['SEXTUPOLE']: sextupole,
    CLASS_CODES['OCTUPOLE']: octupole,
    CLASS_CODES['DECAPOLE']: decapole,
    CLASS_CODES['MULTIPOLE']: multipole,
    CLASS_CODES['HKICKER']: hkicker,
    CLASS_CODES['VKICKER']: vkicker,
    CLASS_CODES['DEGRADER']: degrader,
}
