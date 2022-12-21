# flake8: noqa

import georges_core
import numpy as np
import pytest
from georges_core import ureg

import georges.manzoni
from georges.manzoni import Beam, Input
from georges.manzoni.elements import Quadrupole
from georges.manzoni.maps import (
    compute_mad_combined_dipole_matrix,
    compute_mad_combined_dipole_tensor,
    compute_mad_quadrupole_matrix,
    compute_mad_quadrupole_tensor,
    tmsect,
    track_madx_quadrupole,
)


@pytest.mark.parametrize(
    "h, k1, k2, length, energy",
    [
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
        # Combined function with quadrupolar and sextupolar terms
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
    ],
)
def test_mad8_combined_dipole(h, k1, k2, length, energy):
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


@pytest.mark.parametrize(
    "k1, length, energy, x, px, y, py, dpp",
    [
        # Pure quadrupole
        (1.0, 1.0, 230, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 230, 0.01, 0.0, 0.0, 0.0, 0.0),
        (1.0, 10.0, 230, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.0, 10.0, 230, 0.01, 0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 2300, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 2300, 0.01, 0.0, 0.0, 0.0, 0.0),
        (1.0, 10.0, 2300, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.0, 10.0, 2300, 0.01, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 1.0, 230, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 1.0, 230, 0.01, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 10.0, 230, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 10.0, 230, 0.01, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 1.0, 2300, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 1.0, 2300, 0.01, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 10.0, 2300, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 10.0, 2300, 0.01, 0.0, 0.0, 0.0, 0.0),
    ],
)
def test_mad8_quadrupole(k1, length, energy, x, px, y, py, dpp):
    kin = georges_core.Kinematics(energy * ureg.MeV)

    pt = Beam.compute_pt(
        dpp=dpp,
        beta=kin.beta,
    )

    # MAD-X
    obs_madx = georges.manzoni.BeamObserver()
    input_madx = Input(
        sequence=[
            Quadrupole(
                NAME="Q1",
                L=length * ureg.meter,
                K1=k1 / ureg.meter**2,
                integrator=georges.manzoni.integrators.MadXIntegrator,
            ),
        ],
    )
    beam_madx = Beam(
        kinematics=kin,
        distribution=np.array([[x, px, y, py, dpp, pt]]),
    )
    input_madx.track(
        beam=beam_madx,
        observers=obs_madx,
    )
    res_madx = obs_madx.to_df().at["Q1", "BEAM_OUT"]

    # MAD8
    obs_mad8 = georges.manzoni.BeamObserver()
    input_mad8 = Input(
        sequence=[
            Quadrupole(
                NAME="Q1",
                L=length * ureg.meter,
                K1=k1 / ureg.meter**2,
                integrator=georges.manzoni.integrators.Mad8SecondOrderTaylorIntegrator,
            ),
        ],
    )
    beam_mad8 = Beam(
        kinematics=kin,
        distribution=np.array([[x, px, y, py, 0, pt]]),
    )
    input_mad8.track(
        beam=beam_mad8,
        observers=obs_mad8,
    )
    res_mad8 = obs_mad8.to_df().at["Q1", "BEAM_OUT"]

    assert np.all(np.isclose(res_madx[0][0:4], res_mad8[0][0:4]))
