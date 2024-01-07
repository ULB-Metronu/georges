import matplotlib.pyplot as plt

import georges
from georges import ureg as _ureg
from georges import vis
from georges.manzoni import Input, observers
from georges.manzoni.beam import MadXBeam


def test_manzoni_tracking():
    d1 = georges.Element.Drift(
        NAME="D1",
        L=2.145 * _ureg.m,
        APERTYPE="CIRCULAR",
        APERTURE=[3.645 * _ureg.cm, 3.645 * _ureg.cm],
    )

    d8 = georges.Element.Drift(
        NAME="D8",
        L=0.403 * _ureg.m,
        APERTYPE="CIRCULAR",
        APERTURE=[3.645 * _ureg.cm, 3.645 * _ureg.cm],
    )

    B2G2 = georges.Element.SBend(
        NAME="B2G2",
        L=1.492 * _ureg.m,
        ANGLE=-57 * _ureg.degrees,
        K1=0 * _ureg.m**-2,
        HGAP=0.0 * _ureg.m,
        APERTYPE="RECTANGULAR",
        APERTURE=[5 * _ureg.cm, 2.8 * _ureg.cm],
    )

    d9 = georges.Element.Drift(
        NAME="D9",
        L=0.684 * _ureg.m,
        APERTYPE="CIRCULAR",
        APERTURE=[3.645 * _ureg.cm, 3.645 * _ureg.cm],
    )

    B3G2 = georges.Element.SBend(
        NAME="B3G2",
        L=1.492 * _ureg.m,
        ANGLE=-90 * _ureg.degrees,
        K1=0 * _ureg.m**-2,
        APERTYPE="RECTANGULAR",
        APERTURE=[5 * _ureg.cm, 2.8 * _ureg.cm],
        HGAP=0.0315 * _ureg.m,
        E1=-0.4346 * _ureg.radians,
        E2=-0.2889 * _ureg.radians,
        FINT=0.5,
        FINTX=0.5,
    )

    d10 = georges.Element.Drift(
        NAME="D10",
        L=3.5 * _ureg.m,
        APERTYPE="CIRCULAR",
        APERTURE=[3.645 * _ureg.cm, 3.645 * _ureg.cm],
    )

    sequence = georges.PlacementSequence(name="Sequence")
    sequence.place(d1, at_entry=0 * _ureg.m)
    sequence.place_after_last(d8)
    sequence.place_after_last(B2G2)
    sequence.place_after_last(d9)
    sequence.place_after_last(B3G2)
    sequence.place_after_last(d10)

    kin = georges.Kinematics(250 * _ureg.MeV, particle=georges.particles.Proton, kinetic=True)
    sequence.metadata.kinematics = kin

    beam = MadXBeam(
        kinematics=kin,
        distribution=georges.Distribution.from_twiss_parameters(
            n=100,
            x=2.5 * _ureg.mm,
            y=2.5 * _ureg.mm,
            emitx=7 * _ureg.mm * _ureg.mradians,
            emity=7 * _ureg.mm * _ureg.mradians,
        ).distribution.values,
    )

    mi = Input.from_sequence(sequence=sequence)
    mi.freeze()
    beam_observer_std = mi.track(beam=beam, observers=observers.SigmaObserver())
    beam_observer_beam = mi.track(beam=beam, observers=observers.BeamObserver(with_input_beams=True))
    beam_observer_losses = mi.track(beam=beam, observers=observers.LossesObserver())

    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(111)
    manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
    manzoni_plot.plot_beamline(sequence.df, with_cartouche=True, print_label=True, with_aperture=True)
    manzoni_plot.tracking(beam_observer_std, plane="both")

    assert beam_observer_std is not None
    assert beam_observer_beam is not None
    assert beam_observer_losses is not None
