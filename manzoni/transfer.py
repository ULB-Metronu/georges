import numpy as np
from .constants import *


def rotation(e):
    angle = e[INDEX_ANGLE]
    return np.array(
        [
            [np.cos(angle), 0, -np.sin(angle), 0, 0],
            [0, np.cos(angle), 0, -np.sin(angle), 0],
            [np.sin(angle), 0, np.cos(angle), 0, 0],
            [0, np.sin(angle), 0, np.cos(angle), 0],
            [0, 0, 0, 0, 1],
        ]
    )


def drift(e):
    length = e[INDEX_L]
    return np.array(
        [
            [1, length, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, length, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1]
        ]
    )


def degrader(e):
    NN = 1e5
    np.random_intel.multivariate_normal(
        [1.0, 1.0], np.array(
            [
                [1, 0],
                [0, 1],
            ]),
        int(NN))
    np.random_intel.multivariate_normal(
        [1.0, 1.0], np.array(
            [
                [1, 0],
                [0, 1],
            ]),
        int(NN))
    np.random_intel.normal(0, 1, int(NN))


def sbend(e):
    # http://laacg.lanl.gov/laacg/services/traceman.pdf
    theta = e[INDEX_ANGLE]
    length = e[INDEX_L]
    s = np.sin(theta)
    c = np.cos(theta)
    e1 = e[INDEX_E1]
    e2 = e[INDEX_E2]
    k1 = (-1.0/(length/theta))*np.tan(e1)
    k2 = (-1.0/(length/theta))*np.tan(e2)
    m_e1 = np.array(
        [
            [1, 0, 0, 0, 0],
            [-k1, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, k1, 1, 0],
            [0, 0, 0, 0, 1]
        ]
    )
    m_b = np.array(
        [
            [c, (length / theta) * s, 0, 0, (length/theta)*(1-c)],
            [-(theta / length) * s, c, 0, 0, s],
            [0, 0, 1, length, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1]
        ]
    )
    m_e2 = np.array(
        [
            [1, 0, 0, 0, 0],
            [-k2, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, k2, 1, 0],
            [0, 0, 0, 0, 1]
        ]
    )
    return m_e2 @ m_b @ m_e1


def quadrupole(e):
    length = e[INDEX_L]
    k = e[INDEX_K1]
    if k > 0:
        k = np.sqrt(k)
        kl = k * length
        s = np.sin(kl)
        c = np.cos(kl)
        sh = np.sinh(kl)
        ch = np.cosh(kl)
        return np.array(
            [
                [c, (1 / k) * s, 0, 0, 0],
                [-k * s, c, 0, 0, 0],
                [0, 0, ch, (1 / k) * sh, 0],
                [0, 0, k * sh, ch, 0],
                [0, 0, 0, 0, 1]
            ])
    else:
        k *= -1
        k = np.sqrt(k)
        kl = k * length
        s = np.sin(kl)
        c = np.cos(kl)
        sh = np.sinh(kl)
        ch = np.cosh(kl)
        return np.array(
            [
                [ch, (1 / k) * sh, 0, 0, 0],
                [k * sh, ch, 0, 0, 0],
                [0, 0, c, (1 / k) * s, 0],
                [0, 0, -k * s, c, 0],
                [0, 0, 0, 0, 1]
            ])


transfer = list()
transfer.insert(CLASS_CODE_DRIFT, drift)
transfer.insert(CLASS_CODE_SBEND, sbend)
transfer.insert(CLASS_CODE_QUADRUPOLE, quadrupole)
transfer.insert(CLASS_CODE_NONE, None)
transfer.insert(CLASS_CODE_DEGRADER, degrader)
transfer.insert(CLASS_CODE_ROTATION, rotation)
