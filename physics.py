import numpy as np

PROTON_MASS = 0.938272


def momentum_to_energy(p):
    """Return E [GeV/c^2] from P [GeV/c] (proton)."""
    return 1000*(np.sqrt(p**2+PROTON_MASS**2)-PROTON_MASS)


def momentum_to_brho(p):
    """Return BRHO [T.m] from P [GeV/c] (proton)."""
    return 3.33564 * p


def energy_to_momentum(ekin):
    """Return P [GeV] from E [GeV/c^2] (proton)."""
    E = PROTON_MASS + ekin/1000
    return np.sqrt(E**2-PROTON_MASS**2)


def energy_to_beta(ekin):
    """Return beta relativistic from E [GeV/c^2] (proton)."""
    gamma = (PROTON_MASS + ekin) / PROTON_MASS
    return np.sqrt((gamma ** 2 - 1) / gamma ** 2)


def compute_emittance(x, px):
    """Return the beam emittance [mm.mrad] from x [mm] and px [mrad]."""
    return x * px


def compute_beta(x, px):
    """Return the Twiss beta [m] from beam shape x [mm] and px [mrad]."""
    return x / px


def range_to_energy(r):
    """Return the range [g/cm^2] from the kinetic energy [MeV]."""
    a = 0.00169; b = -0.00490; c = 0.56137; d = 3.46405
    return np.exp(
        a * np.log(r)**3 + b * np.log(r)**2 + c * np.log(r) + d
    )


def energy_to_range(e):
    """Return the kinetic energy [MeV] from the range [g/cm^2]."""
    b = 0.008539; c = 0.5271; d = 3.4917
    return np.exp((-c + np.sqrt(c**2 - 4 * b * (d - np.log(e))))/(2*b))


def compute_ess_transmission(beam_sigma, slits, dispersion):
    """Compute the transmission as a function of the momentum offset (in %) from a simple analytical model."""
    n_steps = 10000
    dx = 3.0/n_steps
    sigma = beam_sigma/2.8
    slits_at = slits/dispersion
    error = np.arange(-1.5, 1.5, dx)
    slits = np.zeros(n_steps)
    slits[np.where((error < slits_at) & (error > -slits_at))] = 1.0
    beam = np.exp(-(error/sigma)**2/2)
    return np.roll(np.convolve(slits, beam, mode="same"), -1)/np.trapz(beam)