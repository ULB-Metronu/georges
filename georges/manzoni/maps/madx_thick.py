from typing import Tuple
from numpy import sin, cos, sinh, cosh, sqrt
from numba import njit, prange


@njit(fastmath=True)
def _apply_tilt_rotation(x: float,
                         xp: float,
                         y: float,
                         yp: float,
                         ct: float,
                         st: float,
                         sign: int) -> Tuple[float, float, float, float]:
    """

    Args:
        x:
        xp:
        y:
        yp:
        ct:
        st:
        sign:

    Returns:

    """
    x_ = ct * x + sign * st * y
    y_ = ct * y - sign * st * x
    px_ = ct * xp + sign * st * yp
    py_ = ct * yp - sign * st * xp
    return x_, px_, y_, py_


@njit(parallel=True, fastmath=True)
def track_madx_srotation(b1, b2, element_parameters: list, global_parameters: list):
    """

    Args:
        b1:
        b2:
        element_parameters:
        global_parameters:

    Returns:

    """
    tilt: float = element_parameters[0]
    if tilt < 1e-8:
        b2 = b1.copy()
    else:
        st = sin(tilt)
        ct = cos(tilt)
        for i in prange(b1.shape[0]):
            x = b1[i, 0]
            px = b1[i, 1]
            y = b1[i, 2]
            py = b1[i, 3]
            x_, px_, y_, py_ = _apply_tilt_rotation(x, px, y, py, ct, st, 1)
            b2[i, 0] = x_
            b2[i, 1] = px_
            b2[i, 2] = y_
            b2[i, 3] = py_
            b2[i, 4] = b1[i, 4]
            b2[i, 5] = b1[i, 5]

    return b1, b2


@njit(parallel=True, fastmath=True)
def track_madx_drift(b1, b2, element_parameters: list, global_parameters: list):
    """
    Track through a drift. This methods follows directly from the method implemented in MAD-X.

    .. note:: From MAD-X: `ttdrf` in `trrun.f90`

    Args:
        b1:
        b2:
        element_parameters:
        global_parameters:

    Returns:

    """
    length: float = element_parameters[0]
    beta: float = global_parameters[0]
    time_of_flight: bool = b1.shape[0] == 7
    for i in prange(b1.shape[0]):
        px = b1[i, 1]
        py = b1[i, 3]
        pt = b1[i, 5]

        lpz = length / sqrt(1.0 + 2.0 * pt / beta + pt**2.0 - px**2.0 - py**2.0)

        b2[i, 0] = b1[i, 0] + lpz * px
        b2[i, 1] = b1[i, 1]
        b2[i, 2] = b1[i, 2] + lpz * py
        b2[i, 3] = b1[i, 3]
        b2[i, 4] = b1[i, 4]
        b2[i, 5] = b1[i, 5]
        if time_of_flight:
            b2[i, 6] = b1[i, 6] + (length - (1.0 + beta * pt) * lpz) / beta

    return b1, b2


@njit(parallel=True, fastmath=True)
def track_madx_drift_paraxial(b1, b2, element_parameters: list, global_parameters: list):
    """
    Track through a drift using the paraxial approximation. Not used by MAD-X.

    Args:
        b1:
        b2:
        element_parameters:
        global_parameters:

    Returns:

    """
    length: float = element_parameters[0]
    beta: float = global_parameters[0]
    time_of_flight: bool = b1.shape[0] == 7
    for i in prange(b1.shape[0]):
        b2[i, 0] = b1[i, 0] + length * b1[i, 1]  # X
        b2[i, 1] = b1[i, 1]  # PX
        b2[i, 2] = b1[i, 2] + length * b1[i, 3]  # Y
        b2[i, 3] = b1[i, 3]  # PY
        b2[i, 4] = b1[i, 4]  # DPP
        b2[i, 5] = b1[i, 5]  # PT
        if time_of_flight:
            b2[i, 6] = b1[i, 6] + (length - (1.0 + beta * b1[i, 5]) * length) / beta

    return b1, b2


@njit(parallel=True, fastmath=True)
def track_madx_quadrupole(b1, b2, element_parameters: list, global_parameters: list):
    """
    Track through a (thick) quadrupole. This methods follows directly from the method implemented in MAD-X.

    The Hamiltonian is

    H = (1/2) K1 x^2 + (1/2) px^2/(delta + 1)

    .. note:: From MAD-X: `ttcfd` in `trrun.f90`

    Args:
        b1:
        b2:
        element_parameters:
        global_parameters:

    Returns:

    """
    length = element_parameters[0]
    k1 = element_parameters[1]
    tilt = element_parameters[2]
    st: float
    ct: float

    if k1 == 0:
        return track_madx_drift(b1, b2, element_parameters, global_parameters)

    if tilt != 0.0:
        st = sin(tilt)
        ct = cos(tilt)

    for i in prange(b1.shape[0]):
        delta_plus_1 = b1[i, 4] + 1
        x = b1[i, 0]
        xp = b1[i, 1] / delta_plus_1  # This is the key point to remember
        y = b1[i, 2]
        yp = b1[i, 3] / delta_plus_1  # This is the key point to remember

        if tilt != 0.0:
            x, xp, y, yp = _apply_tilt_rotation(x, xp, y, yp, ct, st, 1)

        k1_ = k1 / delta_plus_1  # This is the key point to remember
        if k1_ > 0:
            kl = sqrt(k1_) * length
            sx = sin(kl) / sqrt(k1_)
            cx = cos(kl)
            sy = sinh(kl) / sqrt(k1_)
            cy = cosh(kl)
        else:
            kl = sqrt(-k1_) * length
            sx = sinh(kl) / sqrt(-k1_)
            cx = cosh(kl)
            sy = sin(kl) / sqrt(-k1_)
            cy = cos(kl)

        x_ = cx * x + sx * xp
        xp_ = (-k1_ * sx * x + cx * xp) * delta_plus_1
        y_ = cy * y + sy * yp
        yp_ = (k1_ * sy * y + cy * yp) * delta_plus_1

        if tilt != 0.0:
            x_, xp_, y_, yp_ = _apply_tilt_rotation(x_, xp_, y_, yp_, ct, st, -1)

        b2[i, 0] = x_
        b2[i, 1] = xp_
        b2[i, 2] = y_
        b2[i, 3] = yp_
        b2[i, 4] = b1[i, 4]
        b2[i, 5] = b1[i, 5]

    return b1, b2


@njit(parallel=True, fastmath=True)
def track_madx_bend(b1, b2, element_parameters: list, global_parameters: list):
    """
    Track through a (thick) combined function bend. This methods follows directly from the method implemented in MAD-X.

    .. note:: From MAD-X: `ttcfd` in `trrun.f90`

    Args:
        b1:
        b2:
        element_parameters:
        global_parameters:

    Returns:

    """
    length: float = element_parameters[0]
    angle: float = element_parameters[1]
    k1: float = element_parameters[2]
    # k2 = element_parameters[3]
    tilt: float = element_parameters[4]
    h: float = element_parameters[5]
    k0: float = element_parameters[6]
    entrance_fringe_x: float = element_parameters[7]
    entrance_fringe_y: float = element_parameters[8]
    exit_fringe_x: float = element_parameters[9]
    exit_fringe_y: float = element_parameters[10]

    if angle == 0:
        return track_madx_quadrupole(b1,
                                     b2,
                                     [
                                         length, k1, 0.0, tilt
                                     ],
                                     global_parameters)

    if tilt != 0.0:
        st: float = sin(tilt)
        ct: float = cos(tilt)

    for i in prange(b1.shape[0]):
        delta_plus_1: float = b1[i, 4] + 1.0
        x: float = b1[i, 0]
        xp: float = b1[i, 1]
        y: float = b1[i, 2]
        yp: float = b1[i, 3]

        # Apply magnet rotation
        if tilt != 0.0:
            xp, xp, y, yp = _apply_tilt_rotation(x, xp, y, yp, ct, st, 1)

        # Apply entrance fringe field
        xp += entrance_fringe_x * x
        yp += entrance_fringe_y * y

        # Body of the magnet
        k0_ = k0 / delta_plus_1
        k1_ = k1 / delta_plus_1
        # k2_ = k2 / delta_plus_1
        kx = k0_ * h + k1_
        ky = -k1_

        if kx > 0:
            klx = sqrt(kx) * length
            sx = sin(klx) / sqrt(kx)
            cx = cos(klx)
        elif kx < 0:
            klx = sqrt(-kx) * length
            sx = sinh(klx) / sqrt(-kx)
            cx = cosh(klx)
        else:
            sx = length
            cx = 1

        if ky > 0:
            kly = sqrt(ky) * length
            sy = sin(kly) / sqrt(ky)
            cy = cos(kly)
        elif ky < 0:
            kly = sqrt(-ky) * length
            sy = sinh(kly) / sqrt(-ky)
            cy = cosh(kly)
        else:
            sy = length
            cy = 1

        xp /= delta_plus_1
        yp /= delta_plus_1

        x_: float = cx * x + sx * xp
        xp_: float = ((-kx * x - k0_ + h) * sx + cx * xp) * delta_plus_1
        y_: float = cy * y + sy * yp
        yp_: float = (-ky * sy * y + cy * yp) * delta_plus_1

        if kx != 0.0:
            x_ = x_ + (k0_ - h) * (cx - 1.0) / kx
        else:
            x_ = x_ - (k0_ - h) * 0.5 * length**2

        # Apply exit fringe field
        xp_ += exit_fringe_x * x_
        yp_ += exit_fringe_y * y_

        # Apply magnet rotation
        if tilt != 0.0:
            x_, xp_, y_, yp_ = _apply_tilt_rotation(x_, xp_, y_, yp_, ct, st, -1)

        b2[i, 0] = x_
        b2[i, 1] = xp_
        b2[i, 2] = y_
        b2[i, 3] = yp_
        b2[i, 4] = b1[i, 4]
        b2[i, 5] = b1[i, 5]

    return b1, b2


@njit(parallel=True, fastmath=True)
def track_madx_dipedge(b1, b2, element_parameters: list, global_parameters: list):
    fringe_x: float = element_parameters[0]
    fringe_y: float = element_parameters[1]

    for i in prange(b1.shape[0]):
        x: float = b1[i, 0]
        xp: float = b1[i, 1]
        y: float = b1[i, 2]
        yp: float = b1[i, 3]

        # Apply entrance fringe field
        xp += fringe_x * x
        yp += fringe_y * y

        b2[i, 0] = x
        b2[i, 1] = xp
        b2[i, 2] = y
        b2[i, 3] = yp
        b2[i, 4] = b1[i, 4]
        b2[i, 5] = b1[i, 5]

    return b1, b2


@njit(parallel=True, fastmath=True)
def track_madx_kicker(b1, b2, element_parameters: list, global_parameters: list):
    hkick: float = element_parameters[1]
    vkick: float = element_parameters[2]

    # Half drift
    track_madx_drift(b1, b2, element_parameters, global_parameters)

    # Kick
    for i in prange(b1.shape[0]):
        b1[i, 0] = b2[i, 0]
        b1[i, 1] = b2[i, 1] + hkick
        b1[i, 2] = b2[i, 2]
        b1[i, 3] = b2[i, 3] + vkick
        b1[i, 4] = b2[i, 4]
        b1[i, 5] = b2[i, 5]

    # Half drift
    track_madx_drift(b1, b2, element_parameters, global_parameters)

    return b1, b2
