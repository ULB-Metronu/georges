import numpy as np


def get_range_from_energy(material, energy, **kwargs):
    csda = kwargs.get("csda", False)
    projected = kwargs.get("projected", not csda)
    if projected and not csda:
        return np.exp(materials_db['projected_ranges'][material](np.log(energy))) / density[material]
    elif csda and not projected:
        return np.exp(materials_db['csda_ranges'][material](np.log(energy))) / density[material]
    else:
        raise Exception("'projected' or 'csda' arguments are mutually exclusive and one must be defined.")


def get_energy_from_range(material, r, **kwargs):
    csda = kwargs.get("csda", False)
    projected = kwargs.get("projected", not csda)
    if projected and not csda:
        return \
        np.exp(materials_db['projected_ranges'][material].solve(np.log(r * density[material]), extrapolate=False))[0]
    elif csda and not projected:
        return np.exp(materials_db['csda_ranges'][material].solve(np.log(r * density[material]), extrapolate=False))[0]
    else:
        raise Exception("'projected' or 'csda' arguments are mutually exclusive and one must be defined.")


def residual_energy(material, thickness, k_in):
    return get_energy_from_range(material, residual_range(material, thickness, k_in))


def residual_range(material, thickness, k_in):
    return get_range_from_energy(material, k_in) - thickness


def necessary_thickness(material, k_out, k_in, **kwargs):
    return get_range_from_energy(material, k_out, **kwargs) - get_range_from_energy(material, k_in, **kwargs)