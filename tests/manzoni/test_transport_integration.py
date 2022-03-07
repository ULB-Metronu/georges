import pytest
import numpy as np
import georges.manzoni
import georges.manzoni.integrators
from georges.manzoni import *
from georges.manzoni.maps import compute_transport_combined_dipole_matrix, \
    compute_transport_combined_dipole_tensor, compute_transport_quadrupole_matrix, compute_transport_quadrupole_tensor
import georges_core
from georges_core import ureg


@pytest.mark.parametrize(
    "k1, k2, length, energy, x, px, y, py, dpp", [
        # Combined function with K1 (K2 = 0)
        (1.0, 0.0, 1.0, 230, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 0.0, 1.0, 230, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.0, 0.0, 10.0, 230, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 0.0, 10.0, 230, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.0, 0.0, 1.0, 2300, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 0.0, 1.0, 2300, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.0, 0.0, 10.0, 2300, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-1.0, 0.0, 10.0, 2300, 0.0, 0.0, 0.0, 0.0, 0.0),
    ])
def test_transport_combined_dipole( k1, k2, length, energy, x, px, y, py, dpp):
    kin = georges_core.Kinematics(energy * ureg.MeV)

    pt = Beam.compute_pt(dpp=dpp,
                         beta=kin.beta)

    # MAD-X
    obs_madx = georges.manzoni.BeamObserver()
    input_madx = Input(sequence=[SBend(NAME='D1',
                                       L=length * ureg.meter,
                                       ANGLE=30 * ureg.rad,
                                       K1=k1 / ureg.m ** 2,
                                       integrator=georges.manzoni.integrators.MadXIntegrator)
                                 ])
    beam_madx = Beam(kinematics=kin,
                     distribution=np.array([[x, px, y, py, dpp, pt]]))
    input_madx.track(beam=beam_madx,
                     observers=obs_madx)
    res_madx = obs_madx.to_df().at['D1', 'BEAM_OUT']

    # Transport
    obs_transport = georges.manzoni.BeamObserver()
    input_transport = Input(sequence=[SBend(NAME='D1',
                                            L=length * ureg.meter,
                                            ANGLE=30 * ureg.rad,
                                            K1=k1 / ureg.m ** 2,
                                            integrator=georges.manzoni.integrators.TransportSecondOrderTaylorIntegrator)
                                      ])
    beam_transport = Beam(kinematics=kin,
                          distribution=np.array([[x, px, y, py, 0, dpp]]))
    input_transport.track(beam=beam_transport,
                          observers=obs_transport)
    res_transport = obs_transport.to_df().at['D1', 'BEAM_OUT']

    assert np.all(np.isclose(res_madx[0][0:4], res_transport[0][0:4]))


@pytest.mark.parametrize(
    "k1, length, energy, x, px, y, py, dpp", [
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
    ])
def test_transport_quadrupole(k1, length, energy, x, px, y, py, dpp):
    kin = georges_core.Kinematics(energy * ureg.MeV)

    pt = Beam.compute_pt(dpp=dpp,
                         beta=kin.beta)

    # MAD-X
    obs_madx = georges.manzoni.BeamObserver()
    input_madx = Input(sequence=[Quadrupole(NAME='Q1',
                                            L=length * ureg.meter,
                                            K1=k1 / ureg.meter ** 2,
                                            integrator=georges.manzoni.integrators.MadXIntegrator)
                                 ])
    beam_madx = Beam(kinematics=kin,
                     distribution=np.array([[x, px, y, py, dpp, pt]]))
    input_madx.track(beam=beam_madx,
                     observers=obs_madx)
    res_madx = obs_madx.to_df().at['Q1', 'BEAM_OUT']

    # Transport
    obs_transport = georges.manzoni.BeamObserver()
    input_transport = Input(sequence=[Quadrupole(NAME='Q1',
                                                 L=length * ureg.meter,
                                                 K1=k1 / ureg.meter ** 2,
                                                 integrator=georges.manzoni.integrators.TransportSecondOrderTaylorIntegrator)
                                      ])
    beam_transport = Beam(kinematics=kin,
                          distribution=np.array([[x, px, y, py, 0, dpp]]))
    input_transport.track(beam=beam_transport,
                          observers=obs_transport)
    res_transport = obs_transport.to_df().at['Q1', 'BEAM_OUT']

    assert np.all(np.isclose(res_madx[0][0:4], res_transport[0][0:4]))
