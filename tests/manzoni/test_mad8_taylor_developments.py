import pytest
import numpy as np
from georges.manzoni.maps import compute_mad_combined_dipole_matrix, compute_mad_combined_dipole_tensor
from georges.manzoni.maps import tmsect
import georges_core
from georges_core import ureg
from numpy import *

def __compute_mad_combined_dipole_tensor(element_parameters: list, global_parameters: list) -> np.ndarray:
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    K1: float = element_parameters[2]
    K2: float = element_parameters[3]
    beta: float = global_parameters[0]
    gamma = 1 / sqrt(1 - beta ** 2)
    h = alpha / L
    kx2 = h ** 2 + K1
    ky2 = -K1
    T = np.zeros((6, 6, 6))

    # Setup basic notations - Horizontal plane
    if h ** 2 + K1 == 0:
        cx = 1
        sx = L
        dx = L ** 2 / 2
        fx = L ** 3 / 6
        c2x = 1
        s2x = L
        d2x = L ** 2 / 2
        f2x = L ** 3 / 6
    elif h ** 2 + K1 > 0:
        print(f"h**2 + K1 > 0: {h ** 2 + K1}")
        cx = cos(sqrt(h ** 2 + K1) * L)
        sx = sin(sqrt(h ** 2 + K1) * L) / sqrt(h ** 2 + K1)
        dx = (1 - cos(sqrt(h ** 2 + K1) * L)) / (h ** 2 + K1)
        fx = (L - sin(sqrt(h ** 2 + K1) * L) / sqrt(h ** 2 + K1)) / (h ** 2 + K1)
        c2x = cos(2 * sqrt(h ** 2 + K1) * L)
        s2x = sin(2 * sqrt(h ** 2 + K1) * L) / (2 * sqrt(h ** 2 + K1))
        d2x = (1 - cos(2 * sqrt(h ** 2 + K1) * L)) / (4 * (h ** 2 + K1))
        f2x = (L - sin(2 * sqrt(h ** 2 + K1) * L) / (2 * sqrt(h ** 2 + K1))) / (4 * (h ** 2 + K1))
    else:
        cx = cosh(sqrt(-h ** 2 - K1) * L)
        sx = sinh(sqrt(-h ** 2 - K1) * L) / sqrt(-h ** 2 - K1)
        dx = (1 - cosh(sqrt(-h ** 2 - K1) * L)) / (h ** 2 + K1)
        fx = (L - sinh(sqrt(-h ** 2 - K1) * L) / sqrt(-h ** 2 - K1)) / (h ** 2 + K1)
        c2x = cosh(2 * sqrt(-h ** 2 - K1) * L)
        s2x = sinh(2 * sqrt(-h ** 2 - K1) * L) / (2 * sqrt(-h ** 2 - K1))
        d2x = (1 - cosh(2 * sqrt(-h ** 2 - K1) * L)) / (4 * (h ** 2 + K1))
        f2x = (L - sinh(2 * sqrt(-h ** 2 - K1) * L) / (2 * sqrt(-h ** 2 - K1))) / (4 * (h ** 2 + K1))

    # Setup basic notations - Vertical plane
    if K1 == 0:
        print(f"K1 == 0: {K1}")
        cy = 1
        sy = L
        dy = L ** 2 / 2
        fy = L ** 3 / 6
        c2y = 1
        s2y = 1 * L
        d2y = L ** 2 / 2
        f2y = L ** 3 / 6
    elif K1 > 0:
        cy = cosh(sqrt(K1) * L)
        sy = sinh(sqrt(K1) * L) / sqrt(K1)
        dy = -((1 - cosh(sqrt(K1) * L)) / K1)
        fy = -((L - sinh(sqrt(K1) * L) / sqrt(K1)) / K1)
        c2y = cosh(2 * sqrt(K1) * L)
        s2y = sinh(2 * sqrt(K1) * L) / (2 * sqrt(K1))
        d2y = -(1 - cosh(2 * sqrt(K1) * L)) / (4 * K1)
        f2y = -(L - sinh(2 * sqrt(K1) * L) / (2 * sqrt(K1))) / (4 * K1)
    else:
        cy = cos(sqrt(-K1) * L)
        sy = sin(sqrt(-K1) * L) / sqrt(-K1)
        dy = -((1 - cos(sqrt(-K1) * L)) / K1)
        fy = -((L - sin(sqrt(-K1) * L) / sqrt(-K1)) / K1)
        c2y = cos(2 * sqrt(-K1) * L)
        s2y = sin(2 * sqrt(-K1) * L) / (2 * sqrt(-K1))
        d2y = -(1 - cos(2 * sqrt(-K1) * L)) / (4 * K1)
        f2y = -(L - sin(2 * sqrt(-K1) * L) / (2 * sqrt(-K1))) / (4 * K1)

    # Special notations
    if K1 == 0 and h == 0:
        cxp2y = 1
        cxm2y = 1
        sxp2y = L
        sxm2y = L
        dxp2y = L ** 2 / 2
        dxm2y = L ** 2 / 2
        fxp2y = L ** 3 / 6
        fxm2y = L ** 3 / 6
    elif kx2 >= 0.0 and ky2 >= 0:
        print(f"not K1 == 0 and h == 0: K1={K1}, h={h}")
        cxp2y = cos(((2 * sqrt(-K1) + sqrt(h ** 2 + K1)) * L) / 2)
        cxm2y = cos(((-2 * sqrt(-K1) + sqrt(h ** 2 + K1)) * L) / 2)
        sxp2y = (2 * sin(((2 * sqrt(-K1) + sqrt(h ** 2 + K1)) * L) / 2)) / (2 * sqrt(-K1) + sqrt(h ** 2 + K1))
        sxm2y = (2 * sin(((-2 * sqrt(-K1) + sqrt(h ** 2 + K1)) * L) / 2)) / (-2 * sqrt(-K1) + sqrt(h ** 2 + K1))
        dxp2y = (4 * (1 - cos(((2 * sqrt(-K1) + sqrt(h ** 2 + K1)) * L) / 2))) / (
                    2 * sqrt(-K1) + sqrt(h ** 2 + K1)) ** 2
        dxm2y = (4 * (1 - cos(((-2 * sqrt(-K1) + sqrt(h ** 2 + K1)) * L) / 2))) / (
                    -2 * sqrt(-K1) + sqrt(h ** 2 + K1)) ** 2
        fxp2y = (4 * (L - (2 * sin(((2 * sqrt(-K1) + sqrt(h ** 2 + K1)) * L) / 2)) / (
                    2 * sqrt(-K1) + sqrt(h ** 2 + K1)))) / (2 * sqrt(-K1) + sqrt(h ** 2 + K1)) ** 2
        fxm2y = (4 * (L - (2 * sin(((-2 * sqrt(-K1) + sqrt(h ** 2 + K1)) * L) / 2)) / (
                    -2 * sqrt(-K1) + sqrt(h ** 2 + K1)))) / (-2 * sqrt(-K1) + sqrt(h ** 2 + K1)) ** 2

    # Compute the basic integrals - Integrals J1, J2 and J3
    if kx2 * L ** 2 < 1e-2:
        j1 = L ** 3 / 6 - (kx2 * L ** 5) / 120 + (kx2 ** 2 * L ** 7) / 5040
        j2 = L ** 5 / 20 - (kx2 * L ** 7) / 168 + (kx2 ** 2 * L ** 9) / 2880
        j3 = L ** 7 / 56 - (kx2 * L ** 9) / 288 + (7 * kx2 ** 2 * L ** 11) / 21120
    else:
        print(f"kx2 * L**2 > 1e-2: kx2={kx2}")
        j1 = (L - sx) / kx2
        j2 = (3 * L - 4 * sx + cx * sx) / (2 * kx2 ** 2)
        j3 = (15 * L - 22 * sx + 9 * cx * sx - 2 * cx ** 2 * sx) / (6 * kx2 ** 3)

    # Compute the basic integrals - Integrals JC, JS, JD, JF
    jc = 0
    js = 0
    jd = 0
    jf = 0
    print(f"kx2={kx2}")  # 1
    print(f"ky2={ky2}")  # 0
    if max(kx2, 4 * ky2) < 1e-2 or (kx2 - 4 * ky2) < 1e-6:
        print("1")
        jc = L ** 2 / 2 + ((4 * K1 - kx2) * L ** 4) / 24 + ((16 * K1 ** 2 - 4 * K1 * kx2 + kx2 ** 2) * L ** 6) / 720
        js = L ** 3 / 6 + ((4 * K1 - kx2) * L ** 5) / 120 + ((16 * K1 ** 2 - 4 * K1 * kx2 + kx2 ** 2) * L ** 7) / 5040
        jd = L ** 4 / 24 + ((4 * K1 - kx2) * L ** 6) / 720 + ((16 * K1 ** 2 - 4 * K1 * kx2 + kx2 ** 2) * L ** 8) / 40320
        jf = L ** 5 / 120 + ((4 * K1 - kx2) * L ** 7) / 5040 + (
                    (16 * K1 ** 2 - 4 * K1 * kx2 + kx2 ** 2) * L ** 9) / 362880
    elif kx2 <= 0 or ky2 <= 0:
        print("2")
        jc = (c2y - cx) / (kx2 - 4 * ky2)
        js = (s2y - sx) / (kx2 - 4 * ky2)
        jd = (d2y - dx) / (kx2 - 4 * ky2)
        # zd = (dyy - dx) / dd
        # xksq - y2ksq
        jf = (f2y - fx) / (kx2 - 4 * ky2)
    elif kx2 > 4 * ky2 > 0:
        print("3")
        kx = sqrt(kx2)
        ky = sqrt(ky2)
        jc = (sxm2y * sxp2y) / 2
        js = -(cxp2y * sxm2y) / (4 * kx * ky) + cxm2y * sxp2y
        jd = (d2y - (sxm2y * sxp2y) / 2) / kx2
        jf = (f2y + (cxp2y * sxm2y) / (4 * kx * ky) - cxm2y * sxp2y) / kx2
    elif 4 * ky2 > kx2 > 0:
        print("4")
        kx = sqrt(kx2)
        ky = sqrt(ky2)
        jc = (sxm2y * sxp2y) / 2
        js = -(cxp2y * sxm2y) / (4 * kx * ky) + cxm2y * sxp2y
        jd = (dx - (sxm2y * sxp2y) / 2) / ky2
        jf = (fx + (cxp2y * sxm2y) / (4 * kx * ky) - cxm2y * sxp2y) / ky2

    # Definition of the tensor elements
    T[0, 0, 0] = -(h * kx2 * sx ** 2) / 2 - ((2 * h * K1 + K2) * (dx + sx ** 2)) / 6
    T[0, 0, 1] = (cx * h * sx) / 2 - (dx * (2 * h * K1 + K2) * sx) / 6
    T[0, 1, 1] = (cx * dx * h) / 2 - (dx ** 2 * (2 * h * K1 + K2)) / 6
    T[0, 0, 5] = (K1 * L * sx) / (4 * beta) + (h ** 2 * sx ** 2) / (2 * beta) - (
                h * (2 * h * K1 + K2) * (-dx ** 2 + 3 * j1 * sx)) / (12 * beta)
    T[0, 1, 5] = -(cx * L + sx) / (4 * beta) + (h ** 2 * (cx * j1 + dx * sx)) / (4 * beta) - (
                h * (2 * h * K1 + K2) * (-2 * cx * j2 + dx ** 2 * sx)) / (12 * beta)
    T[0, 5, 5] = (h ** 3 * j1 * sx) / (2 * beta ** 2) - (h * L * sx) / (2 * beta ** 2) - (
                h ** 2 * (2 * h * K1 + K2) * (dx ** 3 - 2 * j2 * sx)) / (6 * beta ** 2) - (dx * h) / (
                             2 * beta ** 2 * gamma ** 2)
    T[0, 2, 2] = jd * K1 * K2 + (dx * (h * K1 + K2)) / 2
    T[0, 2, 3] = (js * K2) / 2
    T[0, 3, 3] = -(dx * h) / 2 + jd * K2
    T[1, 0, 0] = -((1 + 2 * cx) * (2 * h * K1 + K2) * sx) / 6
    T[1, 0, 1] = -((1 + 2 * cx) * dx * (2 * h * K1 + K2)) / 6
    T[1, 1, 1] = -(h * sx) / 2 - (dx * (2 * h * K1 + K2) * sx) / 3
    T[1, 0, 5] = -(K1 * (-(cx * L) + sx)) / (4 * beta) - (h * (2 * h * K1 + K2) * (3 * cx * j1 + dx * sx)) / (12 * beta)
    T[1, 1, 5] = (K1 * L * sx) / (4 * beta) - (h * (2 * h * K1 + K2) * (dx ** 2 + 3 * j1 * sx)) / (12 * beta)
    T[1, 5, 5] = -(h * K1 * (cx * j1 - dx * sx)) / (2 * beta ** 2) - (
                h ** 2 * (2 * h * K1 + K2) * (-2 * cx * j2 + dx ** 2 * sx)) / (6 * beta ** 2) - (h * sx) / (
                             2 * beta ** 2 * gamma ** 2)
    T[1, 2, 2] = js * K1 * K2 + ((h * K1 + K2) * sx) / 2
    T[1, 2, 3] = (jc * K2) / 2
    T[1, 3, 3] = js * K2 - (h * sx) / 2
    T[2, 0, 2] = (h * K1 * sx * sy) / 2 + (K2 * (cy * jc - 2 * js * K1 * sy)) / 2
    T[2, 0, 3] = (cy * h * sx) / 2 + (K2 * (-2 * cy * js + jc * sy)) / 2
    T[2, 1, 2] = (dx * h * K1 * sy) / 2 + (K2 * (cy * js - 2 * jd * K1 * sy)) / 2
    T[2, 1, 3] = (cy * dx * h) / 2 + (K2 * (-2 * cy * jd + js * sy)) / 2
    T[2, 2, 5] = (h ** 2 * j1 * K1 * sy) / (2 * beta) - (K1 * L * sy) / (4 * beta) + (
                h * K2 * (cy * jd - 2 * jf * K1 * sy)) / (2 * beta)
    T[2, 3, 5] = (cy * h ** 2 * j1) / (2 * beta) - (cy * L + sy) / (4 * beta) + (h * K2 * (-2 * cy * jf + jd * sy)) / (
                2 * beta)
    T[3, 0, 2] = (cy * (h * K1 + K2) * sx) / 2 + (K1 * K2 * (2 * cy * js - jc * sy)) / 2
    T[3, 0, 3] = ((h * K1 + K2) * sx * sy) / 2 + (K2 * (-(cy * jc) + 2 * js * K1 * sy)) / 2
    T[3, 1, 2] = (cy * dx * (h * K1 + K2)) / 2 + (K1 * K2 * (2 * cy * jd - js * sy)) / 2
    T[3, 1, 3] = (dx * (h * K1 + K2) * sy) / 2 + (K2 * (-(cy * js) + 2 * jd * K1 * sy)) / 2
    T[3, 2, 5] = (cy * h * j1 * (h * K1 + K2)) / (2 * beta) + (K1 * (-(cy * L) + sy)) / (4 * beta) + (
                h * K1 * K2 * (2 * cy * jf - jd * sy)) / (2 * beta)
    T[3, 3, 5] = (h * j1 * (h * K1 + K2) * sy) / (2 * beta) - (K1 * L * sy) / (4 * beta) + (
                h * K2 * (-(cy * jd) + 2 * jf * K1 * sy)) / (2 * beta)
    T[4, 0, 0] = -(K1 * (L - cx * sx)) / (4 * beta) + (h * (2 * h * K1 + K2) * (3 * j1 + dx * sx)) / (12 * beta)
    T[4, 0, 1] = (dx ** 2 * h * (2 * h * K1 + K2)) / (12 * beta) + (K1 * sx ** 2) / (4 * beta)
    T[4, 1, 1] = (h * j2 * (2 * h * K1 + K2)) / (6 * beta) - sx / (2 * beta) - (K1 * (j1 - dx * sx)) / (4 * beta)
    T[4, 0, 5] = ((1 + cx) * h * j1 * K1) / (4 * beta ** 2) + (h ** 2 * (3 * dx * j1 - 4 * j2) * (2 * h * K1 + K2)) / (
                12 * beta ** 2) + (h * sx) / (2 * beta ** 2 * gamma ** 2)
    T[4, 1, 5] = (h * j1 * K1 * sx) / (4 * beta ** 2) + (h ** 2 * (2 * h * K1 + K2) * (dx ** 3 - 2 * j2 * sx)) / (
                12 * beta ** 2) + (dx * h) / (2 * beta ** 2 * gamma ** 2)
    T[4, 5, 5] = (h ** 3 * (-2 * dx * j2 + 3 * j3) * (2 * h * K1 + K2)) / (6 * beta ** 3) + (
                h ** 2 * K1 * (-((1 + 2 * cx) * j2) + dx ** 2 * sx)) / (6 * beta ** 3) + (3 * (h ** 2 * j1 - L)) / (
                             2 * beta ** 3 * gamma ** 2)
    T[4, 2, 2] = -((h * jf * K1 * K2) / beta) - (h * j1 * (h * K1 + K2)) / (2 * beta) + (K1 * (L - cy * sy)) / (
                4 * beta)
    T[4, 2, 3] = -(h * jd * K2) / (2 * beta) - (K1 * sy ** 2) / (4 * beta)
    T[4, 3, 3] = (h ** 2 * j1) / (2 * beta) - (h * jf * K2) / beta - (L + cy * sy) / (4 * beta)

    return T


@pytest.mark.parametrize(
    "h, k1, k2, L, energy", [
        (0.0, 0.0, 0.0, 1.0, 230),
        (1.0, 0.0, 0.0, 1.0, 230),
        (1.0, 1.0, 1.0, 1.0, 230),
    ])
def test_tensor_validity(h, k1, k2, L, energy):
    kin = georges_core.Kinematics(energy * ureg.MeV)

    # MAD-X
    r_madx, t_madx = tmsect(fsec=True, el=L, h=h, sk1=k1, sk2=k2, dh=0, beta=kin.beta, gamma=kin.gamma)

    # MAD8
    element_parameters = [L, L*h, k1, k2]
    global_parameters = [kin.beta]

    r_mad8 = compute_mad_combined_dipole_matrix(element_parameters, global_parameters)
    t_mad8 = compute_mad_combined_dipole_tensor(element_parameters, global_parameters)

    assert np.all(np.isclose(r_madx, r_mad8))
    print(np.isclose(t_madx, t_mad8))
    # print(t_madx)
    # print(t_mad8)
    print(t_madx[2, 2, 5])
    print(t_mad8[2, 2, 5])
    assert np.all(np.isclose(t_madx[0], t_mad8[0]))
    assert np.all(np.isclose(t_madx[1], t_mad8[1]))
    assert np.all(np.isclose(t_madx[2], t_mad8[2]))
    assert np.all(np.isclose(t_madx[3], t_mad8[3]))
    assert np.all(np.isclose(t_madx[4], t_mad8[4]))
    assert np.all(np.isclose(t_madx[5], t_mad8[5]))
