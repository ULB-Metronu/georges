from scipy.integrate import quad
from ..physics import energy_to_pv
from .stopping import residual_energy


def fermi_eyges_integrals(u, initial_energy, thickness, material, db, t, n):
    return (thickness-u)**n * t.t(
        energy_to_pv(residual_energy(material, u, initial_energy, db=db)),
        energy_to_pv(initial_energy),
        db=db,
        material=material
    )


def compute_fermi_eyges(material, energy, thickness, db, t, **kwargs):
    a = [
        quad(fermi_eyges_integrals, 0, thickness, args=(energy, thickness, material, db, t, 0))[0],  # Order 0
        quad(fermi_eyges_integrals, 0, thickness, args=(energy, thickness, material, db, t, 1))[0],  # Order 1
        quad(fermi_eyges_integrals, 0, thickness, args=(energy, thickness, material, db, t, 2))[0],  # Order 2
    ]
    b = a[0] * a[2] - a[1]**2  # Emittance

    return {
        'A': a,
        'B': b,
        'E_R': residual_energy(material, thickness, energy, db=db),
        'DPP': 0.0,  # TODO
    }
