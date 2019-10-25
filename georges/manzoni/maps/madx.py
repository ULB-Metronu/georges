from numpy import sin, cos, sinh, cosh, sqrt
from numba import njit, prange


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
    length = element_parameters[0]
    for i in prange(b1.shape[0]):
        beta = global_parameters[0]
        px = b1[i, 1]
        py = b1[i, 3]
        pt = b1[i, 5]

        lpz = length / sqrt(1 + 2 * pt / beta + pt**2 - px**2 - py**2)
        b2[i, 0] = b1[i, 0] + lpz * px
        b2[i, 2] = b1[i, 2] + lpz * py
        b2[i, 4] = b1[i, 4] + (length - (1 + beta * pt) * lpz) / beta

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
    if k1 == 0:
        return track_madx_drift(b1, b2, element_parameters, global_parameters)

    for i in prange(b1.shape[0]):
        delta_plus_1 = b1[i, 4] + 1
        x = b1[i, 0]
        xp = b1[i, 1] / delta_plus_1  #
        y = b1[i, 2]
        yp = b1[i, 3] / delta_plus_1

        k1_ = k1 / delta_plus_1  # This is the key point to remember
        if k1_ > 0:
            kl = sqrt(k1_) * length
            sx = sin(kl) / sqrt(k1_)
            cx = cos(kl)
            sy = sinh(kl)
            cy = cosh(kl)
        else:
            kl = sqrt(-k1_) * length
            sx = sinh(kl) / sqrt(-k1_)
            cx = cosh(kl)
            sy = sin(kl)
            cy = cos(kl)

        x_ = cx * x + sx * xp
        xp_ = (-(sqrt(k1_) * sx) * x + cx * xp) * delta_plus_1
        y_ = cy * y + sy * yp
        yp_ = (sqrt(k1_) * sy * y + cy * yp) * delta_plus_1

        b2[i, 0] = x_
        b2[i, 1] = xp_
        b2[i, 2] = y_
        b2[i, 3] = yp_
        b2[i, 4] = b1[i, 4]
        b2[i, 5] = b1[i, 5]

    return b1, b2
