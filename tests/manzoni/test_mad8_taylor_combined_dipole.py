import pytest
import numpy as np
from georges.manzoni.maps import compute_mad_combined_dipole_matrix, compute_mad_combined_dipole_tensor
from georges.manzoni.maps import tmsect
import georges_core
from georges_core import ureg


@pytest.mark.parametrize(
    "h, k1, k2, length, energy", [
        # Drift like
        (0.0, 0.0, 0.0, 1.0, 230),
        (0.0, 0.0, 0.0, 10.0, 230),
        (0.0, 0.0, 0.0, 1.0, 2300),
        (0.0, 0.0, 0.0, 10.0, 2300),
        # Pure dipole
        (1.0, 0.0, 0.0, 1.0, 230),
        (1.0, 0.0, 0.0, 10.0, 230),
        (1.0, 0.0, 0.0, 1.0, 2300),
        (1.0, 0.0, 0.0, 10.0, 2300),
        # Pure quadrupole
        (0.0, 1.0, 0.0, 1.0, 230),
        (0.0, 1.0, 0.0, 10.0, 230),
        (0.0, 1.0, 0.0, 1.0, 2300),
        (0.0, 1.0, 0.0, 10.0, 2300),
        (0.0, -1.0, 0.0, 1.0, 230),
        (0.0, -1.0, 0.0, 10.0, 230),
        (0.0, -1.0, 0.0, 1.0, 2300),
        (0.0, -1.0, 0.0, 10.0, 2300),
        # Pure sextupole
        (0.0, 0.0, 1.0, 1.0, 230),
        (0.0, 0.0, 1.0, 10.0, 230),
        (0.0, 0.0, 1.0, 1.0, 2300),
        (0.0, 0.0, 1.0, 10.0, 2300),
        (0.0, 0.0, -1.0, 1.0, 230),
        (0.0, 0.0, -1.0, 10.0, 230),
        (0.0, 0.0, -1.0, 1.0, 2300),
        (0.0, 0.0, -1.0, 10.0, 2300),
        # Combined function with K1 (K2 = 0)
        (1.0, 1.0, 0.0, 1.0, 230),
        (1.0, -1.0, 0.0, 1.0, 230),
        (1.0, 1.0, 0.0, 10.0, 230),
        (1.0, -1.0, 0.0, 10.0, 230),
        (1.0, 1.0, 0.0, 1.0, 2300),
        (1.0, -1.0, 0.0, 1.0, 2300),
        (1.0, 1.0, 0.0, 10.0, 2300),
        (1.0, -1.0, 0.0, 10.0, 2300),
        # Combined function with K2 (K1 = 0)
        (1.0, 0.0, 1.0, 1.0, 230),
        (1.0, 0.0, -1.0, 1.0, 230),
        (1.0, 0.0, 1.0, 10.0, 230),
        (1.0, 0.0, -1.0, 10.0, 230),
        (1.0, 0.0, 1.0, 1.0, 2300),
        (1.0, 0.0, -1.0, 1.0, 2300),
        (1.0, 0.0, 1.0, 10.0, 2300),
        (1.0, 0.0, -1.0, 10.0, 2300),
        # Combined funciton with quadrupolar and sextupolar terms
        (1.0, 1.0, 1.0, 1.0, 230),
        (1.0, 1.0, -1.0, 1.0, 230),
        (1.0, -1.0, 1.0, 1.0, 230),
        (1.0, -1.0, -1.0, 1.0, 230),
        (1.0, 1.0, 1.0, 10.0, 230),
        (1.0, 1.0, -1.0, 10.0, 230),
        (1.0, -1.0, 1.0, 10.0, 230),
        (1.0, -1.0, -1.0, 10.0, 230),
        (1.0, 1.0, 1.0, 1.0, 2300),
        (1.0, 1.0, -1.0, 1.0, 2300),
        (1.0, -1.0, 1.0, 1.0, 2300),
        (1.0, -1.0, -1.0, 1.0, 2300),
        (1.0, 1.0, 1.0, 10.0, 2300),
        (1.0, 1.0, -1.0, 10.0, 2300),
        (1.0, -1.0, 1.0, 10.0, 2300),
        (1.0, -1.0, -1.0, 10.0, 2300),
    ])
def test_mad8_taylor_combined_dipole(h, k1, k2, length, energy):
    kin = georges_core.Kinematics(energy * ureg.MeV)

    # MAD-X
    r_madx, t_madx = tmsect(fsec=True, el=length, h=h, sk1=k1, sk2=k2, dh=0, beta=kin.beta, gamma=kin.gamma)

    # MAD8
    element_parameters = [length, length * h, k1, k2]
    global_parameters = [kin.beta]

    r_mad8 = compute_mad_combined_dipole_matrix(element_parameters, global_parameters)
    t_mad8 = compute_mad_combined_dipole_tensor(element_parameters, global_parameters)

    assert np.all(np.isclose(r_madx, r_mad8))
    assert np.all(np.isclose(t_madx, t_mad8))
