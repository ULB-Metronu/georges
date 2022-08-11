import numpy as np
import sympy as sy


def generate_transport_matrix_tensor(h, k1, k2, order, sign):
    """
        h:
        k1: quadrupolar component
        k2: sextupolar component
        order: 1=matrix ; 2=tensor
        sign: sign of h**2 + k1 (-1, 0, 1)
    """

    a, b, c, d = sy.symbols('a b c d')
    L = sy.symbols('L', positive=True, real=True)
    ky, t, tau = sy.symbols('ky t tau', real=True)

    if sign == -1:
        kx2 = sy.symbols('kx2', real=True, negative=True, nonzero=True)
    elif sign == 0:
        kx2 = 0
    elif sign == 1:
        kx2 = sy.symbols('kx2', real=True, positive=True, nonzero=True)
    else:
        assert sign in [-1, 0, 1], f"sign must be -1, 0 or 1, got: {sign}"

    x1 = sy.Function('x1')(t)
    y1 = sy.Function('y1')(t)

    # General solution of EDs :
    xsol = sy.dsolve(sy.Eq(x1.diff(t, t) + kx2 * x1, 0), x1).rhs
    ysol = sy.dsolve(sy.Eq(y1.diff(t, t) - k1 * y1, 0), y1).rhs

    # Limit conditions :
    cnd0x = sy.Eq(xsol.subs(t, 0), a)  # x(0) = a
    cnd1x = sy.Eq(sy.diff(xsol, t).subs(t, 0), b)  # x'(0) = b

    cnd0y = sy.Eq(ysol.subs(t, 0), c)  # x(0) = c
    cnd1y = sy.Eq(sy.diff(ysol, t).subs(t, 0), d)  # x'(0) = d

    # Solve for C1, C2:
    C1, C2 = sy.symbols('C1, C2')
    C1C2_sl = sy.solve([cnd0x, cnd1x], (C1, C2))

    # Solve for C3, C4:

    C3C4_sl = sy.solve([cnd0y, cnd1y], (C1, C2))

    # Substitute back into solution:

    xsol_0 = sy.simplify(xsol.subs(C1C2_sl))
    ysol_0 = sy.simplify(ysol.subs(C3C4_sl))

    Cx = sy.simplify(xsol_0.subs(a, 1).subs(b, 0).as_real_imag()[0])
    Sx = sy.simplify(xsol_0.subs(a, 0).subs(b, 1).as_real_imag()[0])

    Cy = sy.simplify(ysol_0.subs(c, 1).subs(d, 0).as_real_imag()[0])
    Sy = sy.simplify(ysol_0.subs(c, 0).subs(d, 1).as_real_imag()[0])

    # Green functions :
    Gx = Sx * Cx.subs(t, tau) - Cx * Sx.subs(t, tau)
    Gy = Sy * Cy.subs(t, tau) - Cy * Sy.subs(t, tau)

    # First order :
    if order == 1:

        R = sy.tensor.array.MutableDenseNDimArray.zeros(6, 6)
        R[0, 0] = Cx
        R[0, 1] = Sx
        R[0, 5] = sy.integrate(h * Gx, (tau, 0, t))

        R[1, 0] = Cx.diff(t)
        R[1, 1] = Sx.diff(t)
        R[1, 5] = R[0, 5].diff(t)

        R[2, 2] = Cy
        R[2, 3] = Sy

        R[3, 2] = Cy.diff(t)
        R[3, 3] = Sy.diff(t)

        for i in range(0, 6):
            for j in range(0, 6):
                if R[i, j] != 0:
                    R[i, j] = sy.simplify(R[i, j].subs(t, L).as_real_imag()[0])
                    if sign != 0:
                        R[i, j] = sy.simplify(R[i, j].subs(kx2, h ** 2 + k1))

        return R
    # Second order :
    if order == 2:
        # Driving terms :
        F2 = sy.MutableDenseNDimArray.zeros(6, 6, 6)
        dx = sy.integrate(h * Gx, (tau, 0, t))

        F2[0, 0, 0] = -(h ** 3 + k2 + 2 * k1 * h) * Cx ** 2 + sy.diff(h, t) * Cx * sy.diff(Cx, t) + sy.Rational(1,
                                                                                                                2) * h * (
                          sy.diff(Cx, t)) ** 2
        F2[0, 0, 1] = -2 * (h ** 3 + k2 + 2 * k1 * h) * Cx * Sx + sy.diff(h, t) * (
                Cx * sy.diff(Sx, t) + Sx * sy.diff(Cx, t)) + h * (sy.diff(Cx, t) * sy.diff(Sx, t))
        F2[0, 0, 5] = (2 * h ** 2 + k1) * Cx + -2 * (h ** 3 + k2 + 2 * k1 * h) * Cx * dx + sy.diff(h, t) * (
                Cx * sy.diff(dx, t) + dx * sy.diff(Cx, t)) + h * (sy.diff(Cx, t) * sy.diff(dx, t))
        F2[0, 1, 1] = -(h ** 3 + k2 + 2 * k1 * h) * Sx ** 2 + sy.diff(h, t) * Sx * sy.diff(Sx, t) + sy.Rational(1,
                                                                                                                2) * h * (
                          sy.diff(Sx, t)) ** 2
        F2[0, 1, 5] = (2 * h ** 2 + k1) * Sx + 2 * (h ** 3 + k2 + 2 * k1 * h) * Sx * dx + sy.diff(h, t) * (
                Sx * sy.diff(dx, t) + dx * sy.diff(Sx, t)) + h * (sy.diff(Sx, t) * sy.diff(dx, t))
        F2[0, 5, 5] = -h + (2 * h ** 2 + k1) * dx + (h ** 3 + k2 + 2 * k1 * h) * dx ** 2 + sy.diff(h, t) * dx * sy.diff(
            dx, t) + sy.Rational(1, 2) * h * (sy.diff(dx, t)) ** 2
        F2[0, 2, 2] = sy.Rational(1, 2) * (sy.diff(sy.diff(h, t), t) - k1 * h + 2 * k2) * Cy ** 2 + sy.diff(h,
                                                                                                            t) * Cy * sy.diff(
            Cy, t) - sy.Rational(1, 2) * h * sy.diff(Cy, t) ** 2
        F2[0, 2, 3] = (sy.diff(sy.diff(h, t), t) + k1 * h + 2 * k2) * Cy * Sy + sy.diff(h, t) * (
                Cy * sy.diff(Sy, t) + Sy * sy.diff(Cy, t)) - h * sy.diff(Cy, t) * sy.diff(Sy, t)
        F2[0, 3, 3] = sy.Rational(1, 2) * (sy.diff(sy.diff(h, t), t) + k1 * h + 2 * k2) * Sy ** 2 + sy.diff(h,
                                                                                                            t) * Sy * sy.diff(
            Sy, t) - sy.Rational(1, 2) * h * sy.diff(Sy, t) ** 2

        F2[2, 0, 2] = 2 * (k2 + k1 * h) * Cx * Cy + sy.diff(h, t) * (
                Cx * sy.diff(Cy, t) - sy.diff(Cx, t) * Cy) + h * sy.diff(Cx, t) * sy.diff(Cy, t)
        F2[2, 0, 3] = 2 * (k2 + k1 * h) * Cx * Sy + sy.diff(h, t) * (
                Cx * sy.diff(Sy, t) - sy.diff(Cx, t) * Sy) + h * sy.diff(Cx, t) * sy.diff(Sy, t)
        F2[2, 1, 2] = 2 * (k2 + k1 * h) * Sx * Cy + sy.diff(h, t) * (
                Sx * sy.diff(Cy, t) - sy.diff(Sx, t) * Cy) + h * sy.diff(Sx, t) * sy.diff(Cy, t)
        F2[2, 1, 3] = 2 * (k2 + k1 * h) * Sx * Sy + sy.diff(h, t) * (
                Sx * sy.diff(Sy, t) - sy.diff(Sx, t) * Sy) + h * sy.diff(Sx, t) * sy.diff(Sy, t)
        F2[2, 2, 5] = -k1 * Cy + 2 * (k2 + k1 * h) * Cy * dx - sy.diff(h, t) * (
                Cy * sy.diff(dx, t) - sy.diff(Cy, t) * dx) + h * sy.diff(Cy, t) * sy.diff(dx, t)
        F2[2, 3, 5] = -k1 * Sy + 2 * (k2 + k1 * h) * Cy * dx - sy.diff(h, t) * (
                Sy * sy.diff(dx, t) - sy.diff(Sy, t) * dx) + h * sy.diff(Sy, t) * sy.diff(dx, t)

        # Solution for 2nd order coeficients :
        T = sy.MutableDenseNDimArray.zeros(6, 6, 6)

        T[0, 0, 0] = sy.integrate(F2[0, 0, 0].subs(t, tau) * Gx, (tau, 0, t)).as_real_imag()[0]
        T[0, 0, 1] = sy.integrate(F2[0, 0, 1].subs(t, tau) * Gx, (tau, 0, t)).as_real_imag()[0]
        T[0, 0, 5] = sy.integrate(F2[0, 0, 5].subs(t, tau) * Gx, (tau, 0, t)).as_real_imag()[0]
        T[0, 1, 1] = sy.integrate(F2[0, 1, 1].subs(t, tau) * Gx, (tau, 0, t)).as_real_imag()[0]
        T[0, 1, 5] = sy.integrate(F2[0, 1, 5].subs(t, tau) * Gx, (tau, 0, t)).as_real_imag()[0]
        T[0, 5, 5] = sy.integrate(F2[0, 5, 5].subs(t, tau) * Gx, (tau, 0, t)).as_real_imag()[0]
        T[0, 2, 2] = sy.integrate(F2[0, 2, 2].subs(t, tau) * Gx, (tau, 0, t)).as_real_imag()[0]
        T[0, 2, 3] = sy.integrate(F2[0, 2, 3].subs(t, tau) * Gx, (tau, 0, t)).as_real_imag()[0]
        T[0, 3, 3] = sy.integrate(F2[0, 3, 3].subs(t, tau) * Gx, (tau, 0, t)).as_real_imag()[0]

        T[2, 0, 2] = sy.integrate(F2[2, 0, 2].subs(t, tau) * Gy, (tau, 0, t)).as_real_imag()[0]
        T[2, 0, 3] = sy.integrate(F2[2, 0, 3].subs(t, tau) * Gy, (tau, 0, t)).as_real_imag()[0]
        T[2, 1, 2] = sy.integrate(F2[2, 1, 2].subs(t, tau) * Gy, (tau, 0, t)).as_real_imag()[0]
        T[2, 1, 3] = sy.integrate(F2[2, 1, 3].subs(t, tau) * Gy, (tau, 0, t)).as_real_imag()[0]
        T[2, 2, 5] = sy.integrate(F2[2, 2, 5].subs(t, tau) * Gy, (tau, 0, t)).as_real_imag()[0]
        T[2, 3, 5] = sy.integrate(F2[2, 3, 5].subs(t, tau) * Gy, (tau, 0, t)).as_real_imag()[0]

        T[1, 0, 0] = T[0, 0, 0].diff(t).as_real_imag()[0]
        T[1, 0, 1] = T[0, 0, 1].diff(t).as_real_imag()[0]
        T[1, 0, 5] = T[0, 0, 5].diff(t).as_real_imag()[0]
        T[1, 1, 1] = T[0, 1, 1].diff(t).as_real_imag()[0]
        T[1, 1, 5] = T[0, 1, 5].diff(t).as_real_imag()[0]
        T[1, 5, 5] = T[0, 5, 5].diff(t).as_real_imag()[0]
        T[1, 2, 2] = T[0, 2, 2].diff(t).as_real_imag()[0]
        T[1, 2, 3] = T[0, 2, 3].diff(t).as_real_imag()[0]
        T[1, 3, 3] = T[0, 3, 3].diff(t).as_real_imag()[0]

        T[3, 0, 2] = T[2, 0, 2].diff(t).as_real_imag()[0]
        T[3, 0, 3] = T[2, 0, 3].diff(t).as_real_imag()[0]
        T[3, 1, 2] = T[2, 1, 2].diff(t).as_real_imag()[0]
        T[3, 1, 3] = T[2, 1, 3].diff(t).as_real_imag()[0]
        T[3, 2, 5] = T[2, 2, 5].diff(t).as_real_imag()[0]
        T[3, 3, 5] = T[2, 3, 5].diff(t).as_real_imag()[0]

        if sign != 2:
            for i in range(0, 6):
                for j in range(0, 6):
                    for k in range(0, 6):
                        if T[i, j, k] != 0:
                            T[i, j, k] = T[i, j, k].subs(t, L).subs(kx2, h ** 2 + k1)

        return T


def write_matrix(R, file, ntab=3):
    tab = "\t" * ntab
    for i in range(np.shape(R)[0]):
        for j in range(np.shape(R)[1]):
            if R[i, j] != 0:
                file.write(f"{tab}R[{i}, {j}] = {sy.simplify(R[i, j])}\n")


def write_tensor(T, file, ntab=3):
    tab = "\t" * ntab
    for i in range(np.shape(T)[0]):
        for j in range(np.shape(T)[1]):
            for k in range(np.shape(T)[2]):
                if T[i, j, k] != 0:
                    file.write(f"{tab}T[{i}, {j}, {k}] = {sy.simplify(T[i, j, k])}\n")
