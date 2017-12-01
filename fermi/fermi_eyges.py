from scipy.integrate import quad
from ..physics import energy_to_pv
from .stopping import residual_energy


def fermi_eyges_integrals(u, initial_energy, thickness, material, db, T, n):
    return (thickness-u)**n * T.t(
        energy_to_pv(residual_energy(material, u, initial_energy, db=db)),
        energy_to_pv(initial_energy),
        db=db,
        material=material
    )


def compute_fermi_eyges(**kwargs):
    material = kwargs.get('material')
    energy = kwargs.get('energy')
    thickness = kwargs.get('thickness')
    db = kwargs.get('db')
    T = kwargs.get('T')
    A = [
        quad(fermi_eyges_integrals, 0, thickness, args=(energy, thickness, material, db, T, 0))[0],  # Order 0
        quad(fermi_eyges_integrals, 0, thickness, args=(energy, thickness, material, db, T, 1))[0],  # Order 1
        quad(fermi_eyges_integrals, 0, thickness, args=(energy, thickness, material, db, T, 2))[0],  # Order 2
    ]
    B = A[0] * A[2] - A[1]**2

    return {
        'A': A,
        'B': B,
        'E_R': residual_energy(material, thickness, energy, db=db)
    }
