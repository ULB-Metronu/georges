import numpy as np


def get_range_from_energy(material, energy, **kwargs):
    db = kwargs.get('db', None)
    csda = kwargs.get("csda", False)
    projected = kwargs.get("projected", not csda)
    if projected and not csda:
        return np.exp(db.projected_ranges[material](np.log(energy))) / db.density(material)
    elif csda and not projected:
        return np.exp(db.csda_ranges[material](np.log(energy))) / db.density(material)
    else:
        raise Exception("'projected' or 'csda' arguments are mutually exclusive and one must be defined.")


def get_energy_from_range(material, r, **kwargs):
    db = kwargs.get('db', None)
    csda = kwargs.get("csda", False)
    projected = kwargs.get("projected", not csda)
    if projected and not csda:
        return \
        np.exp(db.projected_ranges[material].solve(np.log(r * db.density(material)), extrapolate=False))[0]
    elif csda and not projected:
        return np.exp(db.csda_ranges[material].solve(np.log(r * db.density(material)), extrapolate=False))[0]
    else:
        raise Exception("'projected' or 'csda' arguments are mutually exclusive and one must be defined.")


def residual_energy(material, thickness, k_in, **kwargs):
    return get_energy_from_range(material,
                                 residual_range(material, thickness, k_in, db=kwargs.get('db', None)),
                                 db=kwargs.get('db', None)
                                 )


def residual_range(material, thickness, k_in, **kwargs):
    return get_range_from_energy(material, k_in, db=kwargs.get('db', None)) - thickness


def necessary_thickness(material, k_out, k_in, **kwargs):
    return get_range_from_energy(material, k_out, db=kwargs.get('db', None)) \
           - get_range_from_energy(material, k_in, db=kwargs.get('db', None))
