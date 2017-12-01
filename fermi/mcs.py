import numpy as np


def f_dm(p1v1, pv):
    return 0.5244 + 0.1975 * np.log10(1 - (pv / p1v1) ** 2) + 0.2320 * np.log10(pv) - 0.0098 * np.log10(pv) * np.log10(
        1 - (pv / p1v1) ** 2)


def t_dm(pv, p1v1):
    es = 15.0  # MeV
    rho = 1.7  # Graphite g/cm3
    alpha = 0.0072973525664  # fine structure constant
    N = 6.022140857e23  # Avogadro
    Z = 6  # Carbon
    re = 2.8179e-15  # classical electron radius
    A = 12  # Carbon
    inv_scattering_length = 0.03372681281618887  # rho * alpha * N * re**2 * (Z**2 / A) * (2 * np.log(33219*(A*Z)**(-1/3))-1)

    return f_dm(p1v1, pv) * (es / pv) ** 2 * (1 / inv_scattering_length)



