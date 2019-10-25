from numpy import sin, cos, sinh, cosh, sqrt
from numba import njit, prange


@njit(parallel=True, fastmath=True)
def track_madx_drift(b1, b2, element_parameters: list, global_parameters: list):
    """
    Track through a (thick) quadrupole. This methods follows directly from the method implemented in MAD-X.

    The Hamiltonian is

    H = (1/2) K1 x^2 + (1/2) px^2/(delta + 1)

    Args:
        b1:
        b2:
        element_parameters:
        global_parameters:

    Returns:

    """
    for i in prange(b1.shape[0]):
        pass

    return b1, b2


@njit(parallel=True, fastmath=True)
def track_madx_quadrupole(b1, b2, element_parameters: list, global_parameters: list):
    """
    Track through a (thick) quadrupole. This methods follows directly from the method implemented in MAD-X.

    The Hamiltonian is

    H = (1/2) K1 x^2 + (1/2) px^2/(delta + 1)

    Args:
        b1:
        b2:
        element_parameters:
        global_parameters:

    Returns:

    """
    length = element_parameters[0]
    k1 = element_parameters[1]
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
