import numpy as np


def drift(e):
    length = e[1]
    return np.array(
        [
            [1, length, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, length, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1]
        ]
    )


def sbend(e):
    t = e['ANGLE']
    l = e['LENGTH']
    s = np.sin(t)
    c = np.cos(t)
    return np.array(
        [
            [c, (l / t) * s, 0, 0, 0],
            [-(t / l) * s, c, 0, 0, 0],
            [0, 0, 1, l, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1]
        ]
    )


def quadrupole(e):
    l = e[1]
    k = e[2]
    kl = k * l
    if k > 0:
        s = np.sin(kl)
        c = np.cos(kl)
        sh = np.sinh(kl)
        ch = np.cosh(kl)
        return np.array(
            [
                [c, (1 / k) * s, 0, 0, 0],
                [-k * s, c, 0, 0, 0],
                [0, 0, ch, (1 / k) * sh, 0],
                [0, 0, -k * sh, ch, 0],
                [0, 0, 0, 0, 1]
            ])
    else:
        k *= -1
        s = np.sin(kl)
        c = np.cos(kl)
        sh = np.sinh(kl)
        ch = np.cosh(kl)
        return np.array(
            [
                [ch, (1 / k) * sh, 0, 0, 0],
                [-k * sh, ch, 0, 0, 0],
                [0, 0, c, (1 / k) * s, 0],
                [0, 0, -k * s, c, 0],
                [0, 0, 0, 0, 1]
            ])


def edge():
    pass


transfer = list()
transfer.insert(0, drift)
transfer.insert(1, None)
transfer.insert(2, quadrupole)
