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
    """
    Transfer matrix of a drift.
    :param e: element definition
    :return: a numpy array representing the 5D transfer matrix
    """
    length = e[INDEX_LENGTH]
    return np.array(
        [
            [1, length, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, length, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1]
        ]
    )


def rbend(e, multiply=True):
    """
    Transfer matrix of a rectangular bend (RBEND); a bend with pole-face angle.
    :param e: element definition
    :param multiply: return
    :return: a numpy array representing the 5D transfer matrix
    """
    return bend(e, e[INDEX_ANGLE]/2, e[INDEX_ANGLE]/2, multiply)


def sbend(e, multiply=True):
    """
    Transfer matrix of a sector bend (SBEND); a bend with no pole-face angle.
    :param e: element definition
    :param multiply: return
    :return: a numpy array representing the 5D transfer matrix
    """
    return bend(e, 0, 0, multiply)


def bend(e, e1, e2, multiply=True):
    """
    Transfer matrix of a generic bend. Pole-face angles are optional.
    :param e: element definition
    :param e1: entrance pole-face angle (added to the element pole-face angle)
    :param e2: exit pole-face angle (added to the element pole-face angle)
    :param multiply: return the
    :return: a numpy array representing the 5D transfer matrix
    """
    # Special case of a drift
    if e[INDEX_ANGLE] == 0 and e[INDEX_K1] == 0:
        return drift(e)

    # Special case of a quadrupole
    if e[INDEX_ANGLE] == 0 and e[INDEX_K1] != 0:
        return quadrupole(e)

    # Definition of the main variables
    theta = e[INDEX_ANGLE]
    length = e[INDEX_LENGTH]
    e1 = e[INDEX_E1] + e1
    e2 = e[INDEX_E2] + e2
    k_bend = (theta/length)**2
    k = np.sqrt(np.abs(k_bend + e[INDEX_K1] / e[INDEX_BRHO]))
    k1 = np.sqrt(np.abs(e[INDEX_K1] / e[INDEX_BRHO]))
    kl = k * length
    k1l = k1 * length
    s_dpp = np.sin(theta)
    c_dpp = np.cos(theta)

    # Horizontal plane
    if k_bend + e[INDEX_K1] / e[INDEX_BRHO] > 0:
        s = np.sin(kl)
        c = np.cos(kl)
        m1 = [c, (1 / k) * s, 0, 0, (length / theta) * (1 - c_dpp)]
        m2 = [-k * s, c, 0, 0, s_dpp]
    elif k_bend + e[INDEX_K1] / e[INDEX_BRHO] < 0:
        sh = np.sinh(kl)
        ch = np.cosh(kl)
        m1 = [ch, (1 / k) * sh, 0, 0, (length / theta) * (1 - c_dpp)]
        m2 = [k * sh, ch, 0, 0, s_dpp]
    else:  # k_bend + e[INDEX_K1] == 0
        m1 = [1, length, 0, 0, (length / theta) * (1 - c_dpp)]
        m2 = [0, 1, 0, 0, s_dpp]

    # Vertical plane
    if e[INDEX_K1] < 0:
        s = np.sin(k1l)
        c = np.cos(k1l)
        m3 = [0, 0, c, (1 / k1) * s, 0]
        m4 = [0, 0, -k1 * s, c, 0]
    elif e[INDEX_K1] > 0:
        sh = np.sinh(k1l)
        ch = np.cosh(k1l)
        m3 = [0, 0, ch, (1 / k1) * sh, 0]
        m4 = [0, 0, k1 * sh, ch, 0]
    else:  # k1 == 0
        m3 = [0, 0, 1, length, 0]
        m4 = [0, 0, 0, 1, 0]

    # Construct 5D matrix
    m_b = np.stack([
        m1,
        m2,
        m3,
        m4,
        [0, 0, 0, 0, 1],
    ])

    # Poleface angle
    if e1 == 0 and e2 == 0:
        return m_b
    else:
        h = theta / length
        psi1 = e[INDEX_FINT] * h * e[INDEX_HGAP] * (1/np.cos(e1)) * (1 + (np.sin(e1))**2)
        psi2 = e[INDEX_FINT] * h * e[INDEX_HGAP] * (1/np.cos(e2)) * (1 + (np.sin(e2))**2)
        k1_x = h * np.tan(e1)
        k1_y = h * np.tan(e1-psi1)
        k2_x = h * np.tan(e2)
        k2_y = h * np.tan(e2-psi2)

        m_e1 = np.array(
            [
                [1, 0, 0, 0, 0],
                [k1_x, 1, 0, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, -k1_y, 1, 0],
                [0, 0, 0, 0, 1]
            ]
        )
        m_e2 = np.array(
            [
                [1, 0, 0, 0, 0],
                [k2_x, 1, 0, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, -k2_y, 1, 0],
                [0, 0, 0, 0, 1]
            ]
        )
        if multiply:
            return m_e2 @ m_b @ m_e1
        else:
            return [m_e2, m_b, m_e1]


def quadrupole(e):
    """
    Quadrupole transfer matrix of an element.
    :param e: element definition
    :return: a numpy array representing the 5D transfer matrix
    """
    length = e[INDEX_LENGTH]
    k = e[INDEX_K1] / e[INDEX_BRHO]
    if k > 0:  # Focusing quadrupole
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
    elif k < 0:  # Defocusing quadrupole
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
    else:  # Zero strengths, this is a drift
        return drift(e)


matrices = {
    CLASS_CODES['DRIFT']: drift,
    CLASS_CODES['COLLIMATOR']: drift,
    CLASS_CODES['SBEND']: sbend,
    CLASS_CODES['RBEND']: rbend,
    CLASS_CODES['QUADRUPOLE']: quadrupole,
    CLASS_CODES['ROTATION']: rotation,
}

