import numpy.testing as npt

import georges
from georges import ureg as _ureg
from georges.manzoni import Input
from georges.manzoni.beam import MadXBeam, TransportBeam
from georges.manzoni.integrators import TransportSecondOrderTaylorIntegratorExact
from georges.manzoni.observers import LossesObserver, MeanObserver


def test_setting_parameters():
    B2G2 = georges.Element.SBend(
        NAME="B2G2",
        L=1.492 * _ureg.m,
        ANGLE=-90 * _ureg.degrees,
        K1=0 * _ureg.m**-2,
        APERTYPE="RECTANGULAR",
        APERTURE=[5 * _ureg.cm, 2.8 * _ureg.cm],
        HGAP=0.0315 * _ureg.m,
        E1=-0.4346 * _ureg.radians,
        E2=-0.2889 * _ureg.radians,
        K0=-0.97 * 1.0528125514711102 / _ureg.m,
        FINT=0.5,
        FINTX=0.5,
    )

    sequence = georges.PlacementSequence(name="Sequence")
    sequence.place(B2G2, at_entry=0 * _ureg.m)
    sequence.place_after_last(B2G2)

    mi = Input.from_sequence(sequence=sequence)
    mi.set_parameters("B2G2", parameters={"K0": 1 * _ureg.m**-1})

    assert mi.get_parameters("B2G2", "K0") == 1 * _ureg.m**-1
    assert mi.get_parameters("B2G2", ["L", "HGAP"]) == {"L": 1.492 * _ureg.m, "HGAP": 0.0315 * _ureg.m}
    assert mi.get_parameters("B2G2") == {
        "NAME": "B2G2",
        "AT_ENTRY": 1.492 * _ureg.m,
        "AT_CENTER": 2.238 * _ureg.m,
        "AT_EXIT": 2.984 * _ureg.m,
        "APERTYPE": "RECTANGULAR",
        "APERTURE": [5 * _ureg.cm, 2.8 * _ureg.cm],
        "KINEMATICS": None,
        "ANGLE": -90 * _ureg.degrees,
        "K0": 1 * _ureg.m**-1,
        "K1": 0 * _ureg.m**-2,
        "K2": 0.0 * _ureg.m**-3,
        "L": 1.492 * _ureg.m,
        "E1": -0.4346 * _ureg.radians,
        "E2": -0.2889 * _ureg.radians,
        "TILT": 0.0 * _ureg.radians,
        "HGAP": 0.0315 * _ureg.m,
        "FINT": 0.5,
        "FINTX": 0.5,
    }


def test_integrator_setting():
    B2G2 = georges.Element.SBend(
        NAME="B2G2",
        L=1.492 * _ureg.m,
        ANGLE=-10 * _ureg.degrees,
        K1=0 * _ureg.m**-2,
        APERTYPE="RECTANGULAR",
        APERTURE=[5 * _ureg.cm, 2.8 * _ureg.cm],
        HGAP=0.0315 * _ureg.m,
        E1=-0.4346 * _ureg.radians,
        E2=-0.2889 * _ureg.radians,
    )

    sequence = georges.PlacementSequence(name="Sequence")
    sequence.place(B2G2, at_entry=0 * _ureg.m)

    kin = georges.Kinematics(250 * _ureg.MeV, particle=georges.particles.Proton, kinetic=True)
    sequence.metadata.kinematics = kin

    mi = Input.from_sequence(sequence=sequence)
    mi.freeze()
    mi.set_integrator(integrator=TransportSecondOrderTaylorIntegratorExact)
    beam = TransportBeam(
        kinematics=kin,
        distribution=georges.Distribution.from_twiss_parameters(
            n=1000,
            x=2.5 * _ureg.mm,
            y=2.5 * _ureg.mm,
            emitx=7 * _ureg.mm * _ureg.mradians,
            emity=7 * _ureg.mm * _ureg.mradians,
        ).distribution.values,
    )
    beam_observer_mean = mi.track(beam=beam, observers=MeanObserver())
    assert beam_observer_mean is not None


def test_degrader_efficiency():
    d1 = georges.Element.Drift(
        NAME="D1",
        L=0.3 * _ureg.m,
    )

    c1 = georges.Element.Degrader(
        NAME="DEG",
        MATERIAL=georges.fermi.materials.Beryllium,
        L=10 * _ureg.cm,
        WITH_LOSSES=True,
    )

    d2 = georges.Element.Drift(
        NAME="D2",
        L=0.3 * _ureg.m,
    )

    sequence = georges.PlacementSequence(name="Sequence")

    sequence.place(d1, at_entry=0 * _ureg.m)
    sequence.place_after_last(c1)
    sequence.place_after_last(d2)

    kin = georges.Kinematics(230 * _ureg.MeV, particle=georges.particles.Proton, kinetic=True)
    sequence.metadata.kinematics = kin

    beam = MadXBeam(
        kinematics=kin,
        distribution=georges.Distribution.from_5d_multigaussian_distribution(
            n=100000,
            xrms=0.1 * _ureg.cm,
            yrms=0.7 * _ureg.cm,
            pxrms=0.01,
            pyrms=0.01,
        ).distribution.values,
    )

    mi = Input.from_sequence(sequence=sequence)
    efficiency = mi.compute_efficiency(input_energy=kin.ekin)

    losses_observer = mi.track(beam=beam, observers=LossesObserver()).to_df()
    efficiency_observer = losses_observer.iloc[-1]["PARTICLES_OUT"] / losses_observer.iloc[0]["PARTICLES_IN"]

    npt.assert_allclose(efficiency_observer, efficiency, rtol=0.02)
