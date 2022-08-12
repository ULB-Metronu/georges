import numpy as np
from numba import njit


@njit(cache=True)
def compute_mad_combined_dipole_matrix(element_parameters: list, global_parameters: list) -> np.ndarray:
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    k1: float = element_parameters[2]
    beta: float = global_parameters[0]
    gamma = 1 / np.sqrt(1 - beta ** 2)
    h = alpha / L
    kx2 = h ** 2 + k1
    ky2 = -k1
    R = np.zeros((6, 6))

    # Setup basic notations - Horizontal plane
    if h ** 2 + k1 == 0:
        cx = 1
        sx = L
        dx = L ** 2 / 2
    elif h ** 2 + k1 > 0:
        cx = np.cos(np.sqrt(h ** 2 + k1) * L)
        sx = np.sin(np.sqrt(h ** 2 + k1) * L) / np.sqrt(h ** 2 + k1)
        dx = (1 - np.cos(np.sqrt(h ** 2 + k1) * L)) / (h ** 2 + k1)
    else:
        cx = np.cosh(np.sqrt(-h ** 2 - k1) * L)
        sx = np.sinh(np.sqrt(-h ** 2 - k1) * L) / np.sqrt(-h ** 2 - k1)
        dx = (1 - np.cosh(np.sqrt(-h ** 2 - k1) * L)) / (h ** 2 + k1)

    # Setup basic notations - Vertical plane
    if k1 == 0:
        cy = 1
        sy = L
        dy = L ** 2 / 2
    elif k1 > 0:
        cy = np.cosh(np.sqrt(k1) * L)
        sy = np.sinh(np.sqrt(k1) * L) / np.sqrt(k1)
        dy = -((1 - np.cosh(np.sqrt(k1) * L)) / k1)
    else:
        cy = np.cos(np.sqrt(-k1) * L)
        sy = np.sin(np.sqrt(-k1) * L) / np.sqrt(-k1)
        dy = -((1 - np.cos(np.sqrt(-k1) * L)) / k1)

    # Compute the basic integrals - Integrals J1, J2 and J3
    if kx2 * L ** 2 < 1e-2:
        j1 = L ** 3 / 6 - (kx2 * L ** 5) / 120 + (kx2 ** 2 * L ** 7) / 5040
    else:
        j1 = (L - sx) / kx2

    # Definition of the matrix elements
    R[0, 0] = cx
    R[0, 1] = sx
    R[0, 5] = (dx * h) / beta
    R[1, 0] = -(kx2 * sx)
    R[1, 1] = cx
    R[1, 5] = (h * sx) / beta
    R[2, 2] = cy
    R[2, 3] = sy
    R[3, 2] = -(ky2 * sy)
    R[3, 3] = cy
    R[4, 0] = -((h * sx) / beta)
    R[4, 1] = -((dx * h) / beta)
    R[4, 4] = 1
    R[4, 5] = -((h ** 2 * j1) / beta ** 2) + L / (beta ** 2 * gamma ** 2)
    R[5, 5] = 1

    return R


@njit(cache=True)
def compute_mad_combined_dipole_tensor(element_parameters: list, global_parameters: list) -> np.ndarray:
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    k1: float = element_parameters[2]
    k2: float = element_parameters[3]
    beta: float = global_parameters[0]
    gamma = 1 / np.sqrt(1 - beta ** 2)
    h = alpha / L
    kx2 = h ** 2 + k1
    ky2 = -k1
    T = np.zeros((6, 6, 6))

    # Setup basic notations - Horizontal plane
    if h ** 2 + k1 == 0:
        cx = 1
        sx = L
        dx = L ** 2 / 2
        fx = L ** 3 / 6
        c2x = 1
        s2x = L
        d2x = L ** 2 / 2
        f2x = L ** 3 / 6
    elif h ** 2 + k1 > 0:
        cx = np.cos(np.sqrt(h ** 2 + k1) * L)
        sx = np.sin(np.sqrt(h ** 2 + k1) * L) / np.sqrt(h ** 2 + k1)
        dx = (1 - np.cos(np.sqrt(h ** 2 + k1) * L)) / (h ** 2 + k1)
        fx = (L - np.sin(np.sqrt(h ** 2 + k1) * L) / np.sqrt(h ** 2 + k1)) / (h ** 2 + k1)
        c2x = np.cos(2 * np.sqrt(h ** 2 + k1) * L)
        s2x = np.sin(2 * np.sqrt(h ** 2 + k1) * L) / (2 * np.sqrt(h ** 2 + k1))
        d2x = (1 - np.cos(2 * np.sqrt(h ** 2 + k1) * L)) / (4 * (h ** 2 + k1))
        f2x = (L - np.sin(2 * np.sqrt(h ** 2 + k1) * L) / (2 * np.sqrt(h ** 2 + k1))) / (4 * (h ** 2 + k1))
    else:
        cx = np.cosh(np.sqrt(-h ** 2 - k1) * L)
        sx = np.sinh(np.sqrt(-h ** 2 - k1) * L) / np.sqrt(-h ** 2 - k1)
        dx = (1 - np.cosh(np.sqrt(-h ** 2 - k1) * L)) / (h ** 2 + k1)
        fx = (L - np.sinh(np.sqrt(-h ** 2 - k1) * L) / np.sqrt(-h ** 2 - k1)) / (h ** 2 + k1)
        c2x = np.cosh(2 * np.sqrt(-h ** 2 - k1) * L)
        s2x = np.sinh(2 * np.sqrt(-h ** 2 - k1) * L) / (2 * np.sqrt(-h ** 2 - k1))
        d2x = (1 - np.cosh(2 * np.sqrt(-h ** 2 - k1) * L)) / (4 * (h ** 2 + k1))
        f2x = (L - np.sinh(2 * np.sqrt(-h ** 2 - k1) * L) / (2 * np.sqrt(-h ** 2 - k1))) / (4 * (h ** 2 + k1))

    # Setup basic notations - Vertical plane
    if k1 == 0:
        cy = 1
        sy = L
        dy = L ** 2 / 2
        fy = L ** 3 / 6
        c2y = 1
        s2y = L
        d2y = L ** 2 / 2
        f2y = L ** 3 / 6
    elif k1 > 0:
        cy = np.cosh(np.sqrt(k1) * L)
        sy = np.sinh(np.sqrt(k1) * L) / np.sqrt(k1)
        dy = -((1 - np.cosh(np.sqrt(k1) * L)) / k1)
        fy = -((L - np.sinh(np.sqrt(k1) * L) / np.sqrt(k1)) / k1)
        c2y = np.cosh(2 * np.sqrt(k1) * L)
        s2y = np.sinh(2 * np.sqrt(k1) * L) / (2 * np.sqrt(k1))
        d2y = -(1 - np.cosh(2 * np.sqrt(k1) * L)) / (4 * k1)
        f2y = -(L - np.sinh(2 * np.sqrt(k1) * L) / (2 * np.sqrt(k1))) / (4 * k1)
    else:
        cy = np.cos(np.sqrt(-k1) * L)
        sy = np.sin(np.sqrt(-k1) * L) / np.sqrt(-k1)
        dy = -((1 - np.cos(np.sqrt(-k1) * L)) / k1)
        fy = -((L - np.sin(np.sqrt(-k1) * L) / np.sqrt(-k1)) / k1)
        c2y = np.cos(2 * np.sqrt(-k1) * L)
        s2y = np.sin(2 * np.sqrt(-k1) * L) / (2 * np.sqrt(-k1))
        d2y = -(1 - np.cos(2 * np.sqrt(-k1) * L)) / (4 * k1)
        f2y = -(L - np.sin(2 * np.sqrt(-k1) * L) / (2 * np.sqrt(-k1))) / (4 * k1)

    # Special notations
    if k1 == 0 and h == 0:
        cxp2y = 1
        cxm2y = 1
        sxp2y = L
        sxm2y = L
        dxp2y = L ** 2 / 2
        dxm2y = L ** 2 / 2
        fxp2y = L ** 3 / 6
        fxm2y = L ** 3 / 6
    elif kx2 >= 0.0 and ky2 >= 0:
        cxp2y = np.cos(((2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) * L) / 2)
        cxm2y = np.cos(((-2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) * L) / 2)
        sxp2y = (2 * np.sin(((2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) * L) / 2)) / (2 * np.sqrt(-k1) +
                                                                                     np.sqrt(h ** 2 + k1))
        sxm2y = (2 * np.sin(((-2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) * L) / 2)) / (-2 * np.sqrt(-k1) +
                                                                                      np.sqrt(h ** 2 + k1))
        dxp2y = (4 * (1 - np.cos(((2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) * L) / 2))) / (
                    2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) ** 2
        dxm2y = (4 * (1 - np.cos(((-2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) * L) / 2))) / (
                    -2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) ** 2
        fxp2y = (4 * (L - (2 * np.sin(((2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) * L) / 2)) / (
                    2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)))) / (2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) ** 2
        fxm2y = (4 * (L - (2 * np.sin(((-2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) * L) / 2)) / (
                    -2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)))) / (-2 * np.sqrt(-k1) + np.sqrt(h ** 2 + k1)) ** 2

    # Compute the basic integrals - Integrals J1, J2 and J3
    if abs(kx2 * L ** 2) < 1e-2:
        j1 = L ** 3 / 6 - (kx2 * L ** 5) / 120 + (kx2 ** 2 * L ** 7) / 5040
        j2 = L ** 5 / 20 - (kx2 * L ** 7) / 168 + (kx2 ** 2 * L ** 9) / 2880
        j3 = L ** 7 / 56 - (kx2 * L ** 9) / 288 + (7 * kx2 ** 2 * L ** 11) / 21120
    else:
        j1 = (L - sx) / kx2
        j2 = (3 * L - 4 * sx + cx * sx) / (2 * kx2 ** 2)
        j3 = (15 * L - 22 * sx + 9 * cx * sx - 2 * cx ** 2 * sx) / (6 * kx2 ** 3)

    # Compute the basic integrals - Integrals JC, JS, JD, JF
    if max(kx2, 4 * ky2) < 1e-2 or abs(kx2 - 4 * ky2) < 1e-6:
        jc = L ** 2 / 2 + ((4 * k1 - kx2) * L ** 4) / 24 + ((16 * k1 ** 2 - 4 * k1 * kx2 + kx2 ** 2) * L ** 6) / 720
        js = L ** 3 / 6 + ((4 * k1 - kx2) * L ** 5) / 120 + ((16 * k1 ** 2 - 4 * k1 * kx2 + kx2 ** 2) * L ** 7) / 5040
        jd = L ** 4 / 24 + ((4 * k1 - kx2) * L ** 6) / 720 + ((16 * k1 ** 2 - 4 * k1 * kx2 + kx2 ** 2) * L ** 8) / 40320
        jf = L ** 5 / 120 + ((4 * k1 - kx2) * L ** 7) / 5040 + (
                    (16 * k1 ** 2 - 4 * k1 * kx2 + kx2 ** 2) * L ** 9) / 362880
    elif kx2 <= 0 or ky2 <= 0:
        jc = (c2y - cx) / (kx2 - 4 * ky2)
        js = (s2y - sx) / (kx2 - 4 * ky2)
        jd = (d2y - dx) / (kx2 - 4 * ky2)
        jf = (f2y - fx) / (kx2 - 4 * ky2)
    elif kx2 > 4 * ky2 > 0:
        kx = np.sqrt(kx2)
        ky = np.sqrt(ky2)
        jc = (sxm2y * sxp2y) / 2
        js = -(cxp2y * sxm2y) / (4 * kx * ky) + cxm2y * sxp2y
        jd = (d2y - (sxm2y * sxp2y) / 2) / kx2
        jf = (f2y + (cxp2y * sxm2y) / (4 * kx * ky) - cxm2y * sxp2y) / kx2
    elif 4 * ky2 > kx2 > 0:
        kx = np.sqrt(kx2)
        ky = np.sqrt(ky2)
        jc = (sxm2y * sxp2y) / 2
        js = -(cxp2y * sxm2y) / (4 * kx * ky) + cxm2y * sxp2y
        jd = (dx - (sxm2y * sxp2y) / 2) / ky2
        jf = (fx + (cxp2y * sxm2y) / (4 * kx * ky) - cxm2y * sxp2y) / ky2

    # Definition of the tensor elements
    T[0, 0, 0] = -(h * kx2 * sx ** 2) / 2 - ((2 * h * k1 + k2) * (dx + sx ** 2)) / 6
    T[0, 0, 1] = (cx * h * sx) / 2 - (dx * (2 * h * k1 + k2) * sx) / 6
    T[0, 1, 1] = (cx * dx * h) / 2 - (dx ** 2 * (2 * h * k1 + k2)) / 6
    T[0, 0, 5] = (k1 * L * sx) / (4 * beta) + (h ** 2 * sx ** 2) / (2 * beta) - (
                h * (2 * h * k1 + k2) * (-dx ** 2 + 3 * j1 * sx)) / (12 * beta)
    T[0, 1, 5] = -(cx * L + sx) / (4 * beta) + (h ** 2 * (cx * j1 + dx * sx)) / (4 * beta) - (
                h * (2 * h * k1 + k2) * (-2 * cx * j2 + dx ** 2 * sx)) / (12 * beta)
    T[0, 5, 5] = (h ** 3 * j1 * sx) / (2 * beta ** 2) - (h * L * sx) / (2 * beta ** 2) - (
                h ** 2 * (2 * h * k1 + k2) * (dx ** 3 - 2 * j2 * sx)) / (6 * beta ** 2) - (dx * h) / (
                             2 * beta ** 2 * gamma ** 2)
    T[0, 2, 2] = jd * k1 * k2 + (dx * (h * k1 + k2)) / 2
    T[0, 2, 3] = (js * k2) / 2
    T[0, 3, 3] = -(dx * h) / 2 + jd * k2
    T[1, 0, 0] = -((1 + 2 * cx) * (2 * h * k1 + k2) * sx) / 6
    T[1, 0, 1] = -((1 + 2 * cx) * dx * (2 * h * k1 + k2)) / 6
    T[1, 1, 1] = -(h * sx) / 2 - (dx * (2 * h * k1 + k2) * sx) / 3
    T[1, 0, 5] = -(k1 * (-(cx * L) + sx)) / (4 * beta) - (h * (2 * h * k1 + k2) * (3 * cx * j1 + dx * sx)) / (12 * beta)
    T[1, 1, 5] = (k1 * L * sx) / (4 * beta) - (h * (2 * h * k1 + k2) * (dx ** 2 + 3 * j1 * sx)) / (12 * beta)
    T[1, 5, 5] = -(h * k1 * (cx * j1 - dx * sx)) / (2 * beta ** 2) - (
                h ** 2 * (2 * h * k1 + k2) * (-2 * cx * j2 + dx ** 2 * sx)) / (6 * beta ** 2) - (h * sx) / (
                             2 * beta ** 2 * gamma ** 2)
    T[1, 2, 2] = js * k1 * k2 + ((h * k1 + k2) * sx) / 2
    T[1, 2, 3] = (jc * k2) / 2
    T[1, 3, 3] = js * k2 - (h * sx) / 2
    T[2, 0, 2] = (h * k1 * sx * sy) / 2 + (k2 * (cy * jc - 2 * js * k1 * sy)) / 2
    T[2, 0, 3] = (cy * h * sx) / 2 + (k2 * (-2 * cy * js + jc * sy)) / 2
    T[2, 1, 2] = (dx * h * k1 * sy) / 2 + (k2 * (cy * js - 2 * jd * k1 * sy)) / 2
    T[2, 1, 3] = (cy * dx * h) / 2 + (k2 * (-2 * cy * jd + js * sy)) / 2
    T[2, 2, 5] = (h ** 2 * j1 * k1 * sy) / (2 * beta) - (k1 * L * sy) / (4 * beta) + (
                h * k2 * (cy * jd - 2 * jf * k1 * sy)) / (2 * beta)
    T[2, 3, 5] = (cy * h ** 2 * j1) / (2 * beta) - (cy * L + sy) / (4 * beta) + (h * k2 * (-2 * cy * jf + jd * sy)) / (
                2 * beta)
    T[3, 0, 2] = (cy * (h * k1 + k2) * sx) / 2 + (k1 * k2 * (2 * cy * js - jc * sy)) / 2
    T[3, 0, 3] = ((h * k1 + k2) * sx * sy) / 2 + (k2 * (-(cy * jc) + 2 * js * k1 * sy)) / 2
    T[3, 1, 2] = (cy * dx * (h * k1 + k2)) / 2 + (k1 * k2 * (2 * cy * jd - js * sy)) / 2
    T[3, 1, 3] = (dx * (h * k1 + k2) * sy) / 2 + (k2 * (-(cy * js) + 2 * jd * k1 * sy)) / 2
    T[3, 2, 5] = (cy * h * j1 * (h * k1 + k2)) / (2 * beta) + (k1 * (-(cy * L) + sy)) / (4 * beta) + (
                h * k1 * k2 * (2 * cy * jf - jd * sy)) / (2 * beta)
    T[3, 3, 5] = (h * j1 * (h * k1 + k2) * sy) / (2 * beta) - (k1 * L * sy) / (4 * beta) + (
                h * k2 * (-(cy * jd) + 2 * jf * k1 * sy)) / (2 * beta)
    T[4, 0, 0] = -(k1 * (L - cx * sx)) / (4 * beta) + (h * (2 * h * k1 + k2) * (3 * j1 + dx * sx)) / (12 * beta)
    T[4, 0, 1] = (dx ** 2 * h * (2 * h * k1 + k2)) / (12 * beta) + (k1 * sx ** 2) / (4 * beta)
    T[4, 1, 1] = (h * j2 * (2 * h * k1 + k2)) / (6 * beta) - sx / (2 * beta) - (k1 * (j1 - dx * sx)) / (4 * beta)
    T[4, 0, 5] = ((1 + cx) * h * j1 * k1) / (4 * beta ** 2) + (h ** 2 * (3 * dx * j1 - 4 * j2) * (2 * h * k1 + k2)) / (
                12 * beta ** 2) + (h * sx) / (2 * beta ** 2 * gamma ** 2)
    T[4, 1, 5] = (h * j1 * k1 * sx) / (4 * beta ** 2) + (h ** 2 * (2 * h * k1 + k2) * (dx ** 3 - 2 * j2 * sx)) / (
                12 * beta ** 2) + (dx * h) / (2 * beta ** 2 * gamma ** 2)
    T[4, 5, 5] = (h ** 3 * (-2 * dx * j2 + 3 * j3) * (2 * h * k1 + k2)) / (6 * beta ** 3) + (
                h ** 2 * k1 * (-((1 + 2 * cx) * j2) + dx ** 2 * sx)) / (6 * beta ** 3) + (3 * (h ** 2 * j1 - L)) / (
                             2 * beta ** 3 * gamma ** 2)
    T[4, 2, 2] = -((h * jf * k1 * k2) / beta) - (h * j1 * (h * k1 + k2)) / (2 * beta) + (k1 * (L - cy * sy)) / (
                4 * beta)
    T[4, 2, 3] = -(h * jd * k2) / (2 * beta) - (k1 * sy ** 2) / (4 * beta)
    T[4, 3, 3] = (h ** 2 * j1) / (2 * beta) - (h * jf * k2) / beta - (L + cy * sy) / (4 * beta)

    return T
