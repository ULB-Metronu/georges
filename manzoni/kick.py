import numpy as np
from .constants import *


try:
    import numpy.random_intel as nprandom
except ModuleNotFoundError:
    import numpy.random as nprandom


def sextupole(e, n, **kwargs):
    return np.array(
        [
            0,
            0,
            0,
            0,
            0
        ]
    )


def octupole(e, n, **kwargs):
    return np.array(
        [
            0,
            0,
            0,
            0,
            0
        ]
    )


def decapole(e, n, **kwargs):
    return np.array(
        [
            0,
            0,
            0,
            0,
            0
        ]
    )


def multipole(e, n, **kwargs):
    return np.array(
        [
            0,
            0,
            0,
            0,
            0
        ]
    )


def hkicker(e, n, **kwargs):
    return np.array(
        [
            0,
            e[INDEX_K1],
            0,
            0,
            0
        ]
    )


def vkicker(e, n, **kwargs):
    return np.array(
        [
            0,
            0,
            0,
            e[INDEX_K1],
            0
        ]
    )


def degrader(e, n, **kwargs):
    s11 = kwargs.get('deg', {}).get('A', [0, 0, 0])[0]
    s12 = kwargs.get('deg', {}).get('A', [0, 0, 0])[1]
    s22 = kwargs.get('deg', {}).get('A', [0, 0, 0])[2]
    dpp = kwargs.get('deg', {}).get('DPP', 0)

    return nprandom.multivariate_normal(
        [0.0, 0.0, 0.0, 0.0, 0.0], np.array(
            [
                [s11, s12, 0, 0, 0],
                [s12, s22, 0, 0, 0],
                [0, 0, s11, s12, 0],
                [0, 0, s12, s22, 0],
                [0, 0, 0, 0, dpp]
            ]),
        int(n))


kick = {
    CLASS_CODES['SEXTUPOLE']: sextupole,
    CLASS_CODES['OCTUPOLE']: octupole,
    CLASS_CODES['DECAPOLE']: decapole,
    CLASS_CODES['MULTIPOLE']: multipole,
    CLASS_CODES['HKICKER']: hkicker,
    CLASS_CODES['VKICKER']: vkicker,
    CLASS_CODES['DEGRADER']: degrader,
}
