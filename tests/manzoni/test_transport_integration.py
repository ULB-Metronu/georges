from itertools import product

import numpy as np
import pytest

import georges.manzoni
import georges.manzoni.integrators
from georges import ureg as _ureg
from georges.manzoni import Input
from georges.manzoni.beam import MadXBeam, TransportBeam
from georges.manzoni.integrators import (
    TransportFirstOrderTaylorIntegrator,
    TransportFirstOrderTaylorIntegratorExact,
    TransportSecondOrderTaylorIntegrator,
    TransportSecondOrderTaylorIntegratorExact,
)


@pytest.mark.parametrize(
    "k1, k2, length, angle, e1, e2, energy, x, px, y, py, dpp, integrator",
    product(
        [-1, 0, 1],  # k1
        [-0.1, 0, 0.1],  # k2
        [1, 10],  # length
        [-10, 0.0, 10],  # angle
        [-0.4346, 0.0, 0.4346],  # e1
        [-0.2889, 0.0, 0.2889],  # e2
        [230, 2300],  # energy
        [0.0],  # x
        [0.0],  # px
        [0.0],  # y
        [0.0],  # py
        [0],  # dpp
        [
            TransportFirstOrderTaylorIntegrator,
            TransportFirstOrderTaylorIntegratorExact,
            TransportSecondOrderTaylorIntegrator,
            TransportSecondOrderTaylorIntegratorExact,
        ],  # Transport integrator
    ),
)
def test_transport_combined_dipole(k1, k2, length, angle, e1, e2, energy, x, px, y, py, dpp, integrator):
    kin = georges.Kinematics(energy * _ureg.MeV)
    sbend = georges.Element.SBend(
        NAME="SB1",
        L=length * _ureg.m,
        ANGLE=angle * _ureg.degrees,
        K1=k1 * _ureg.m**-2,
        K2=k2 * _ureg.m**-3,
        E1=e1 * _ureg.radians,
        E2=e2 * _ureg.radians,
    )

    sequence = georges.PlacementSequence(name="Sequence")
    sequence.place(sbend, at_entry=0 * _ureg.m)

    input_manzoni = Input.from_sequence(sequence=sequence)
    input_manzoni.freeze()

    # MADX Tracking
    beam_madx = MadXBeam(
        kinematics=kin,
        distribution=np.array([[x, px, y, py, dpp]]),
    )
    mean_obs_madx = input_manzoni.track(
        beam=beam_madx,
        observers=georges.manzoni.BeamObserver(),
    )
    res_madx = mean_obs_madx.to_df().at["SB1", "BEAM_OUT"]

    # Transport
    beam_transport = TransportBeam(
        kinematics=kin,
        distribution=np.array([[x, px, y, py, dpp]]),
    )

    input_manzoni.set_integrator(integrator=integrator)

    mean_obs_transport = input_manzoni.track(
        beam=beam_transport,
        observers=georges.manzoni.BeamObserver(),
    )
    res_transport = mean_obs_transport.to_df().at["SB1", "BEAM_OUT"]

    assert np.all(np.isclose(res_madx[0][0:4], res_transport[0][0:4]))


@pytest.mark.parametrize(
    "k1, length, energy, x, px, y, py, dpp, integrator",
    product(
        [-1, 0, 1],  # k1
        [1, 10],  # length
        [230, 2300],  # energy
        [-0.001, 0.0, 0.001],  # x
        [-0.001, 0.0, 0.001],  # px
        [-0.001, 0.0, 0.001],  # y
        [-0.001, 0.0, 0.001],  # py
        [0],  # dpp
        [
            TransportFirstOrderTaylorIntegrator,
            TransportFirstOrderTaylorIntegratorExact,
            TransportSecondOrderTaylorIntegrator,
            TransportSecondOrderTaylorIntegratorExact,
        ],  # Transport integrator
    ),
)
def test_transport_quadrupole(k1, length, energy, x, px, y, py, dpp, integrator):
    kin = georges.Kinematics(energy * _ureg.MeV)
    quad = georges.Element.Quadrupole(
        NAME="Q1",
        L=length * _ureg.m,
        K1=k1 * _ureg.m**-2,
    )

    sequence = georges.PlacementSequence(name="Sequence")
    sequence.place(quad, at_entry=0 * _ureg.m)

    input_manzoni = Input.from_sequence(sequence=sequence)
    input_manzoni.freeze()

    # MADX Tracking
    beam_madx = MadXBeam(
        kinematics=kin,
        distribution=np.array([[x, px, y, py, dpp]]),
    )
    mean_obs_madx = input_manzoni.track(
        beam=beam_madx,
        observers=georges.manzoni.BeamObserver(),
    )
    res_madx = mean_obs_madx.to_df().at["Q1", "BEAM_OUT"]

    # Transport
    beam_transport = TransportBeam(
        kinematics=kin,
        distribution=np.array([[x, px, y, py, dpp]]),
    )

    input_manzoni.set_integrator(integrator=integrator)

    mean_obs_transport = input_manzoni.track(
        beam=beam_transport,
        observers=georges.manzoni.BeamObserver(),
    )
    res_transport = mean_obs_transport.to_df().at["Q1", "BEAM_OUT"]

    assert np.all(np.isclose(res_madx[0][0:4], res_transport[0][0:4]))


@pytest.mark.parametrize(
    "k2, length, energy, x, px, y, py, dpp, integrator",
    product(
        [-0.1, 0, 0.1],  # k2
        [1, 10],  # length
        [230, 2300],  # energy
        [-0.001, 0.0, 0.001],  # x
        [-0.001, 0.0, 0.001],  # px
        [-0.001, 0.0, 0.001],  # y
        [-0.001, 0.0, 0.001],  # py
        [0],  # dpp
        [
            TransportFirstOrderTaylorIntegrator,
            TransportFirstOrderTaylorIntegratorExact,
            TransportSecondOrderTaylorIntegrator,
            TransportSecondOrderTaylorIntegratorExact,
        ],  # Transport integrator
    ),
)
def test_transport_sextupole(k2, length, energy, x, px, y, py, dpp, integrator):
    kin = georges.Kinematics(energy * _ureg.MeV)
    sext = georges.Element.Sextupole(
        NAME="S1",
        L=length * _ureg.m,
        K2=k2 * _ureg.m**-3,
    )

    sequence = georges.PlacementSequence(name="Sequence")
    sequence.place(sext, at_entry=0 * _ureg.m)

    input_manzoni = Input.from_sequence(sequence=sequence)
    input_manzoni.freeze()

    # Transport
    beam_transport = TransportBeam(
        kinematics=kin,
        distribution=np.array([[x, px, y, py, dpp]]),
    )

    input_manzoni.set_integrator(integrator=integrator)

    mean_obs_transport = input_manzoni.track(
        beam=beam_transport,
        observers=georges.manzoni.BeamObserver(),
    )
    res_transport = mean_obs_transport.to_df().at["S1", "BEAM_OUT"]

    assert res_transport is not None


@pytest.mark.parametrize(
    "k1, k2, length, energy, x, px, y, py, dpp, integrator",
    product(
        [-0.1, 0, 0.1],  # k2
        [-0.01, 0, 0.01],  # k2
        [1, 10],  # length
        [230, 2300],  # energy
        [-0.001, 0.0, 0.001],  # x
        [-0.001, 0.0, 0.001],  # px
        [-0.001, 0.0, 0.001],  # y
        [-0.001, 0.0, 0.001],  # py
        [0],  # dpp
        [
            TransportFirstOrderTaylorIntegrator,
            TransportFirstOrderTaylorIntegratorExact,
            TransportSecondOrderTaylorIntegrator,
            TransportSecondOrderTaylorIntegratorExact,
        ],  # Transport integrator
    ),
)
def test_transport_multipole(k1, k2, length, energy, x, px, y, py, dpp, integrator):
    kin = georges.Kinematics(energy * _ureg.MeV)
    mult = georges.Element.Multipole(
        NAME="M1",
        L=length * _ureg.m,
        K1=k1 * _ureg.m**-2,
        K2=k2 * _ureg.m**-3,
    )

    sequence = georges.PlacementSequence(name="Sequence")
    sequence.place(mult, at_entry=0 * _ureg.m)

    input_manzoni = Input.from_sequence(sequence=sequence)
    input_manzoni.freeze()

    # Transport
    beam_transport = TransportBeam(
        kinematics=kin,
        distribution=np.array([[x, px, y, py, dpp]]),
    )

    input_manzoni.set_integrator(integrator=integrator)

    mean_obs_transport = input_manzoni.track(
        beam=beam_transport,
        observers=georges.manzoni.BeamObserver(),
    )
    res_transport = mean_obs_transport.to_df().at["M1", "BEAM_OUT"]

    assert res_transport is not None


@pytest.mark.parametrize(
    "length, angle, e1, e2, hgap, k1, r1, r2, energy, x, px, y, py, dpp, integrator",
    product(
        [1],  # length
        [-10, 0, 10],  # angle
        [-5, 0, 5],  # e1
        [-5, 0, 5],  # e2
        [0, 2],  # hgap
        [-1, 0, 1],  # k1
        [1],  # r1
        [1],  # r2
        [230, 2300],  # energy
        [-0.001, 0.0, 0.001],  # x
        [-0.001, 0.0, 0.001],  # px
        [-0.001, 0.0, 0.001],  # y
        [-0.001, 0.0, 0.001],  # py
        [0],  # dpp
        [
            TransportFirstOrderTaylorIntegrator,
            TransportFirstOrderTaylorIntegratorExact,
            TransportSecondOrderTaylorIntegrator,
            TransportSecondOrderTaylorIntegratorExact,
        ],  # Transport integrator
    ),
)
def test_transport_fringe(length, angle, e1, e2, hgap, k1, r1, r2, energy, x, px, y, py, dpp, integrator):
    kin = georges.Kinematics(energy * _ureg.MeV)
    fin = georges.Element.Fringein(
        NAME="FIN",
        L=length * _ureg.mm,
        ANGLE=angle * _ureg.degrees,
        E1=e1 * _ureg.degrees,
        HGAP=hgap * _ureg.cm,
        R1=r1 * _ureg.cm,
        K1=k1 * _ureg.m**-2,
    )
    fout = georges.Element.Fringeout(
        NAME="FOUT",
        L=length * _ureg.mm,
        ANGLE=angle * _ureg.degrees,
        E2=e2 * _ureg.degrees,
        HGAP=hgap * _ureg.cm,
        R2=r2 * _ureg.cm,
        K1=k1 * _ureg.m**-2,
    )

    sbend = georges.Element.SBend(
        NAME="M1",
        L=1 * _ureg.m,
        ANGLE=angle * _ureg.degrees,
        K1=k1 * _ureg.m**-2,
    )

    sequence = georges.PlacementSequence(name="Sequence")
    sequence.place(fin, at_entry=0 * _ureg.m)
    sequence.place_after_last(sbend)
    sequence.place_after_last(fout)

    input_manzoni = Input.from_sequence(sequence=sequence)
    input_manzoni.freeze()

    # Transport
    beam_transport = TransportBeam(
        kinematics=kin,
        distribution=np.array([[x, px, y, py, dpp]]),
    )

    input_manzoni.set_integrator(integrator=integrator)

    mean_obs_transport = input_manzoni.track(
        beam=beam_transport,
        observers=georges.manzoni.BeamObserver(),
    )
    res_transport = mean_obs_transport.to_df().at["M1", "BEAM_OUT"]

    assert res_transport is not None
