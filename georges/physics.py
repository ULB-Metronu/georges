import numpy as _np

PROTON_MASS = 938.2720813
ELECTRON_MASS = 0.5109989461
CHARGED_PION_MASS = 139.57018
NEUTRAL_PION_MASS = 134.9766


def kinematics(units="MeV", **kwargs):
    """Return a dictionary with all kinematics parameters from a given input."""
    r = 0
    e = 0
    pc = 0
    brho = 0
    beta = 0
    gamma = 1
    if len(kwargs) > 1:
        raise Exception("A single keyword argument is expected: range, energy, momentum, beta, gamma).")
    if kwargs.get("range"):
        r = kwargs.get('range')
        e = range_to_energy(r)
        pc = energy_to_momentum(e)
        brho = momentum_to_brho(pc)
        beta = energy_to_beta(e)
        gamma = beta_to_gamma(beta)
    elif kwargs.get('energy'):
        e = kwargs.get('energy')
        r = energy_to_range(e)
        pc = energy_to_momentum(e)
        brho = momentum_to_brho(pc)
        beta = energy_to_beta(e)
        gamma = beta_to_gamma(beta)
    elif kwargs.get('momentum'):
        pc = kwargs.get('momentum')
        e = momentum_to_energy(pc)
        r = energy_to_range(e)
        brho = momentum_to_brho(pc)
        beta = energy_to_beta(e)
        gamma = beta_to_gamma(beta)
    elif kwargs.get('beta'):
        beta = kwargs.get('beta')
        e = beta_to_energy(beta)
        pc = energy_to_momentum(e)
        r = energy_to_range(e)
        brho = momentum_to_brho(pc)
        gamma = beta_to_gamma(beta)
    elif kwargs.get('gamma'):
        gamma = kwargs.get('gamma')
        e = gamma_to_energy(gamma)
        pc = energy_to_momentum(e)
        brho = momentum_to_brho(pc)
        r = energy_to_range(e)
        beta = energy_to_beta(e)

    return {
        'range': r,
        'energy': e,
        'momentum': pc,
        'brho': brho,
        'beta': beta,
        'gamma': gamma,
    }


def momentum_to_energy(p):
    """Return E [MeV/c^2] from P [MeV/c] (proton)."""
    return _np.sqrt(p ** 2 + PROTON_MASS ** 2) - PROTON_MASS


def momentum_to_brho(p):
    """Return BRHO [T.m] from P [MeV/c] (proton)."""
    return 3.33564E-3 * p


def energy_to_brho(e):
    """Return BRHO [T.m] from E [MeV] (proton)."""
    return 3.33564E-3 * energy_to_momentum(e)


def energy_to_momentum(ekin):
    """Return P [MeV/c] from E [MeV/c^2] (proton)."""
    E = PROTON_MASS + ekin
    return _np.sqrt(E ** 2 - PROTON_MASS ** 2)


def energy_to_beta(ekin):
    """Return beta relativistic from E [MeV/c^2] (proton)."""
    gamma = (PROTON_MASS + ekin) / PROTON_MASS
    return _np.sqrt((gamma ** 2 - 1) / gamma ** 2)


def beta_to_gamma(beta):
    """Return gamma relativistic from beta."""
    return 1/(_np.sqrt(1 - beta ** 2))


def gamma_to_energy(gamma):
    """Return relativistic energy from gamma."""
    return gamma * PROTON_MASS - PROTON_MASS


def beta_to_energy(beta):
    """Return relativistic energy from beta."""
    return beta_to_gamma(beta) * PROTON_MASS - PROTON_MASS


def energy_to_pv(energy):
    """Return relativistic factor 'pv' from kinetic energy (MeV)."""
    E = energy + PROTON_MASS
    return (E**2 - PROTON_MASS**2) / E


def range_to_energy(r):
    """Return the kinetic energy [MeV] from the range [g/cm^2]."""
    a = 0.00169; b = -0.00490; c = 0.56137; d = 3.46405
    return _np.exp(
        a * _np.log(r) ** 3 + b * _np.log(r) ** 2 + c * _np.log(r) + d
    )


def energy_to_range(e):
    """Return the range [g/cm^2] from the kinetic energy [MeV]."""
    """IEC60601 energy to range in water"""
    b = 0.008539; c = 0.5271; d = 3.4917
    return _np.exp((-c + _np.sqrt(c ** 2 - 4 * b * (d - _np.log(e)))) / (2 * b))


def compute_ess_transmission(beam_sigma, slits, dispersion):
    """Compute the transmission as a function of the momentum offset (in %) from a simple analytical model."""
    n_steps = 10000
    dx = 3.0/n_steps
    sigma = beam_sigma/2.8
    slits_at = slits/dispersion
    error = _np.arange(-1.5, 1.5, dx)
    slits = _np.zeros(n_steps)
    slits[_np.where((error < slits_at) & (error > -slits_at))] = 1.0
    beam = _np.exp(-(error / sigma) ** 2 / 2)
    return _np.roll(_np.convolve(slits, beam, mode="same"), -1) / _np.trapz(beam)
