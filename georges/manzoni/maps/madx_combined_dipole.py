"""
This code is from MAD-X twiss.f90. It is used to allow comparisons between the MAD-X tensor and
the generated code for MAD8 or Transport style tensor computation.
"""
import numpy as np
from numpy import cos, cosh, sin, sinh, sqrt

zero = 0
one = 1
two = 2
three = 3
four = 4
six = 6
nine = 9
twelve = 12
fifteen = 15
twty = 20
twty2 = 22
twty4 = 24
thty = 30
foty2 = 42
fvty6 = 56
svty2 = 72
httwty = 120
c1 = one
c2 = one / 2
c3 = one / twty4
c4 = one / 720
s1 = one
s2 = one / 6
s3 = one / httwty
s4 = one / 5040
cg0 = one / twty
cg1 = 5 / 840
cg2 = 21 / 60480
ch0 = one / fvty6
ch1 = 14 / 4032
ch2 = 147 / 443520


def tmfoc(el, sk1):
    #   !----------------------------------------------------------------------*
    #   ! Purpose:                                                             *
    #   !   Compute linear focussing functions.                                *
    #   ! Input:                                                               *
    #   !   el        (double)  element length.                                *
    #   !   sk1       (double)  quadrupole strength.                           *
    #   ! Output:                                                              *
    #   !   c         (double)  cosine-like function.             c(k,l)       *
    #   !   s         (double)  sine-like function.               s(k,l)       *
    #   !   d         (double)  dispersion function.              d(k,l)       *
    #   !   f         (double)  integral of dispersion function.  f(k,l)       *
    #   !----------------------------------------------------------------------*
    #   double precision, intent(IN)  :: el, sk1
    #   double precision, intent(OUT) :: c, s, d, f

    #   double precision :: qk, qkl, qkl2

    #   double precision, parameter :: twty=20d0, thty=30d0, foty2=42d0
    qk = sqrt(abs(sk1))
    qkl = qk * el
    qkl2 = sk1 * el**2
    if abs(qkl2) < 1e-2:
        c = one - qkl2 * (one - qkl2 / twelve) / two
        s = (one - qkl2 * (one - qkl2 / twty) / six) * el
        d = (one - qkl2 * (one - qkl2 / thty) / twelve) * el**2 / two
        f = (one - qkl2 * (one - qkl2 / foty2) / twty) * el**3 / six
    else:
        if qkl2 > zero:
            c = cos(qkl)
            s = sin(qkl) / qk
        else:
            c = cosh(qkl)
            s = sinh(qkl) / qk
        d = (one - c) / sk1
        f = (el - s) / sk1
    return c, s, d, f


def tmsect(fsec, el, h, sk1, sk2, dh, beta, gamma):
    # use twissbeamfi, only : beta, gamma, dtbyds
    # use matrices, only: EYE
    # use math_constfi, only :
    # implicit none
    # !----------------------------------------------------------------------*
    # !     Purpose:                                                         *
    # !     TRANSPORT map for a sector dipole without fringe fields.         *
    # !     Input:                                                           *
    # !     fsec      (logical) if true, return second order terms.          *
    # !     el        (double)  element length.                              *
    # !     h         (double)  reference curvature of magnet.               *
    # !     dh        (double)  dipole field error.                          *
    # !     sk1       (double)  quadrupole strength.                         *
    # !     sk2       (double)  sextupole strengh.                           *
    # !     Output:                                                          *
    # !     ek(6)     (double)  kick due to dipole.                          *
    # !     re(6,6)   (double)  transfer matrix.                             *
    # !     te(6,6,6) (double)  second order terms.                          *
    # !----------------------------------------------------------------------*
    # logical, intent(IN) :: fsec
    # double precision :: el, h, dh, sk1, sk2
    # double precision :: ek(6), re(6,6), te(6,6,6)

    # double precision :: bi, bi2, bi2gi2
    # double precision :: cm, cp, cx, cy, cyy, dd, difsq, dm, dp, dx, dyy
    # double precision :: fm, fp, fx, fyy, gx, h2, hx, sm, sp, sumsq, sx, sy, syy
    # double precision :: t1, t116, t126, t166, t2, t216, t226, t266
    # double precision :: t336, t346, t436, t446, t5, t516, t526, t566
    # double precision :: xk, xkl, xklsq, xksq, xs6
    # double precision :: y0, y1, y2, y2klsq, y2ksq, yk, ykl, yklsq, yksq, ys2
    # double precision :: zc, zd, zf, zs

    # !---- Initialize.
    ek = np.zeros(6)
    re = np.eye(6)
    te = np.zeros((6, 6, 6))

    bi = 1 / beta
    bi2 = bi * bi
    bi2gi2 = one / (beta * gamma) ** 2

    # !---- Horizontal.
    xksq = h**2 + sk1
    xk = sqrt(abs(xksq))
    xkl = xk * el
    xklsq = xksq * el**2

    if abs(xklsq) < 1e-2:
        cx = c1 - xklsq * (c2 - xklsq * c3)
        sx = (s1 - xklsq * (s2 - xklsq * s3)) * el
        dx = (c2 - xklsq * (c3 - xklsq * c4)) * el**2
        fx = (s2 - xklsq * (s3 - xklsq * s4)) * el**3
        gx = (cg0 - xklsq * (cg1 - xklsq * cg2)) * el**5
        hx = (ch0 - xklsq * (ch1 - xklsq * ch2)) * el**7
    else:
        if xklsq > 0:
            cx = cos(xkl)
            sx = sin(xkl) / xk
        else:
            cx = cosh(xkl)
            sx = sinh(xkl) / xk

        dx = (one - cx) / xksq
        fx = (el - sx) / xksq
        gx = (three * el - sx * (four - cx)) / (two * xksq**2)
        hx = (fifteen * el - sx * (twty2 - nine * cx + two * cx**2)) / (six * xksq**3)

    re[0, 0] = cx
    re[0, 1] = sx
    re[0, 5] = h * dx * bi
    re[1, 0] = -xksq * sx
    re[1, 1] = cx
    re[1, 5] = h * sx * bi
    re[4, 1] = -re[0, 5]
    re[4, 0] = -re[1, 5]
    re[4, 5] = el * bi2gi2 - h**2 * fx * bi2

    ek[0] = -dh * dx
    ek[1] = -dh * sx
    ek[4] = h * dh * fx * bi + el  # *dtbyds

    # !---- Vertical.
    yksq = -sk1
    yk = sqrt(abs(yksq))
    ykl = yk * el
    yklsq = yksq * el**2

    if abs(yklsq) < 1e-2:
        cy = c1 - yklsq * (c2 - yklsq * c3)
        sy = (s1 - yklsq * (s2 - yklsq * s3)) * el
    elif yklsq > 0:
        cy = cos(ykl)
        sy = sin(ykl) / yk
    else:
        cy = cosh(ykl)
        sy = sinh(ykl) / yk

    re[2, 2] = cy
    re[2, 3] = sy
    re[3, 2] = -yksq * sy
    re[3, 3] = cy

    ek[2] = 0
    ek[3] = 0

    # !---- Second-order terms.
    if fsec:
        # !---- Pure horizontal terms.
        xs6 = (sk2 + two * h * sk1) / six
        ys2 = (sk2 + h * sk1) / two
        h2 = h / two

        t116 = xs6 * (three * sx * fx - dx**2) - h * sx**2
        t126 = xs6 * (sx * dx**2 - two * cx * gx) - h * sx * dx
        t166 = xs6 * (dx**3 - two * sx * gx) - h2 * dx**2
        t216 = xs6 * (three * cx * fx + sx * dx)
        t226 = xs6 * (three * sx * fx + dx**2)
        t266 = xs6 * (sx * dx**2 - two * cx * gx)
        t516 = h * xs6 * (three * dx * fx - four * gx) + (sk1 / two) * (fx + sx * dx)
        t526 = h * xs6 * (dx**3 - two * sx * gx) + (sk1 / two) * dx**2
        t566 = h * xs6 * (three * hx - two * dx * gx) + (sk1 / two) * gx - fx

        t1 = (sk1 / two) * (dx**2 - sx * fx) - dx
        t2 = (sk1 / two) * (el * dx - fx)
        t5 = fx - sk1 * (gx - fx * dx / two)

        te[0, 0, 0] = -xs6 * (sx**2 + dx) - h2 * xksq * sx**2
        te[0, 0, 1] = (-xs6 * dx + h2 * cx) * sx
        te[0, 1, 1] = (-xs6 * dx + h2 * cx) * dx
        te[0, 0, 5] = (-h2 * t116 + (sk1 / four) * el * sx) * bi
        te[0, 1, 5] = (-h2 * t126 + (sk1 / four) * (el * dx - fx) - sx / two) * bi
        te[0, 5, 5] = (-(h**2) * t166 + h * t1) * bi2 - h2 * dx * bi2gi2
        te[1, 0, 0] = -xs6 * (one + two * cx) * sx
        te[1, 0, 1] = -xs6 * (one + two * cx) * dx
        te[1, 1, 1] = -(two * xs6 * dx + h2) * sx
        te[1, 0, 5] = (-h2 * t216 - (sk1 / four) * (sx - el * cx)) * bi
        te[1, 1, 5] = (-h2 * t226 + (sk1 / four) * el * sx) * bi
        te[1, 5, 5] = (-(h**2) * t266 + h * t2) * bi2 - h2 * sx * bi2gi2

        te[4, 0, 0] = (h2 * xs6 * (sx * dx + three * fx) - (sk1 / four) * (el - cx * sx)) * bi
        te[4, 0, 1] = (h2 * xs6 * dx**2 + (sk1 / four) * sx**2) * bi
        te[4, 1, 1] = (h * xs6 * gx - sk1 * (fx - sx * dx) / four - sx / two) * bi
        te[4, 0, 5] = h2 * ((t516 - sk1 * (el * dx - fx) / two) * bi2 + sx * bi2gi2)
        te[4, 1, 5] = h2 * ((t526 - sk1 * (dx**2 - sx * fx) / two) * bi2 + dx * bi2gi2)
        te[4, 5, 5] = (h**2 * (t566 + t5) * bi2 + (three / two) * (h**2 * fx - el) * bi2gi2) * bi

        # !---- Mixed terms.
        y2ksq = four * yksq
        cyy, syy, dyy, fyy = tmfoc(el, y2ksq)
        y2klsq = y2ksq * el**2
        if max(abs(y2klsq), abs(xklsq)) < 1e-2:
            y0 = one
            y1 = xklsq + y2klsq
            y2 = xklsq**2 + xklsq * y2klsq + y2klsq**2
            zc = (y0 - (y1 - y2 / thty) / twelve) * el**2 / two
            zs = (y0 - (y1 - y2 / foty2) / twty) * el**3 / six
            zd = (y0 - (y1 - y2 / fvty6) / thty) * el**4 / twty4
            zf = (y0 - (y1 - y2 / svty2) / foty2) * el**5 / httwty
        elif xksq <= zero or yksq <= zero:
            dd = xksq - y2ksq
            zc = (cyy - cx) / dd
            zs = (syy - sx) / dd
            zd = (dyy - dx) / dd
            zf = (fyy - fx) / dd
        else:
            sumsq = (xk / two + yk) ** 2
            difsq = (xk / two - yk) ** 2
            cp, sp, dp, fp = tmfoc(el, sumsq)
            cm, sm, dm, fm = tmfoc(el, difsq)
            zc = sp * sm / two
            zs = (sp * cm - cp * sm) / (four * xk * yk)
            if xksq > y2ksq:
                zd = (dyy - zc) / xksq
                zf = (fyy - zs) / xksq
            else:
                zd = (dx - zc) / y2ksq
                zf = (fx - zs) / y2ksq

        t336 = sk2 * (cy * zd - two * sk1 * sy * zf) + h * sk1 * fx * sy
        t346 = sk2 * (sy * zd - two * cy * zf) + h * fx * cy
        t436 = two * ys2 * fx * cy - sk2 * sk1 * (sy * zd - two * cy * zf)
        t446 = two * ys2 * fx * sy - sk2 * (cy * zd - two * sk1 * sy * zf)

        te[0, 2, 2] = +sk2 * sk1 * zd + ys2 * dx
        te[0, 2, 3] = +sk2 * zs / two
        te[0, 3, 3] = +sk2 * zd - h2 * dx
        te[1, 2, 2] = +sk2 * sk1 * zs + ys2 * sx
        te[1, 2, 3] = +sk2 * zc / two
        te[1, 3, 3] = +sk2 * zs - h2 * sx
        te[2, 0, 2] = +sk2 * (cy * zc / two - sk1 * sy * zs) + h2 * sk1 * sx * sy
        te[2, 0, 3] = +sk2 * (sy * zc / two - cy * zs) + h2 * sx * cy
        te[2, 1, 2] = +sk2 * (cy * zs / two - sk1 * sy * zd) + h2 * sk1 * dx * sy
        te[2, 1, 3] = +sk2 * (sy * zs / two - cy * zd) + h2 * dx * cy
        te[2, 2, 5] = (h2 * t336 - sk1 * el * sy / four) * bi
        te[2, 3, 5] = (h2 * t346 - (sy + el * cy) / four) * bi
        te[3, 0, 2] = sk2 * sk1 * (cy * zs - sy * zc / two) + ys2 * sx * cy
        te[3, 0, 3] = sk2 * (sk1 * sy * zs - cy * zc / two) + ys2 * sx * sy
        te[3, 1, 2] = sk2 * sk1 * (cy * zd - sy * zs / two) + ys2 * dx * cy
        te[3, 1, 3] = sk2 * (sk1 * sy * zd - cy * zs / two) + ys2 * dx * sy
        te[3, 2, 5] = (h2 * t436 + sk1 * (sy - el * cy) / four) * bi
        te[3, 3, 5] = (h2 * t446 - sk1 * el * sy / four) * bi

        te[4, 2, 2] = (-h * sk2 * sk1 * zf - h * ys2 * fx + sk1 * (el - cy * sy) / four) * bi
        te[4, 2, 3] = (-h * sk2 * zd / two - sk1 * sy**2 / four) * bi
        te[4, 3, 3] = (-h * sk2 * zf + h * h2 * fx - (el + sy * cy) / four) * bi
        # call tmsymm(te)

    return re, te
