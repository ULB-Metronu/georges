import numpy as np
import cpymad.madx
from georges import ureg as _
from georges import manzoni
from georges_core import Kinematics
from georges.manzoni.integrators import MadXIntegrator


def test_sbend():
    k = Kinematics(230 * _.MeV)
    mi = manzoni.Input(sequence=[
        manzoni.SBend('B1', integrator=MadXIntegrator, L=2.0 * _.m, ANGLE=1 * _.radian, K1=1 * _.m ** -2),
    ]).freeze()
    pt = 0.02
    deltap = np.sqrt(pt ** 2 + 2 * pt / k.beta + 1) - 1
    distribution = b = np.array([
        [0.1, 0.1, 0.1, 0.1, deltap, pt]
    ])
    beam = manzoni.Beam(kinematics=k, distribution=distribution)
    observer = manzoni.BeamObserver()
    manzoni.track(beamline=mi, beam=beam, observer=observer)
    r = observer.to_df().iloc[-1]['BEAM_OUT']

    b = np.array([
        [0.1, 0.1, 0.1, 0.1, 0.02]
    ])
    m = cpymad.madx.Madx(stdout=False)
    m.input(f"""
    OPTION, SYNRAD=false;
    OPTION, RBARC=false;
    BEAM, PARTICLE=PROTON, BETA={k.beta};

    SEQ: SEQUENCE, REFER=ENTRY, L=2.0;
      B1: SBEND, AT=0.0, L=2.0, ANGLE=1, K1=1;
      M1: MARKER, AT=2.0;
    ENDSEQUENCE;

    USE, SEQUENCE=SEQ;
    """)

    # #
    # D1: DRIFT, AT=0.0, L=1.0;
    #   Q1: QUADRUPOLE, AT=1.0, L=1.6, K1=-16, TILT=0.2;

    m.command.track(DELTAP=0.0, ONEPASS=True, ONETABLE=True, QUANTUM=False)

    for i in range(b.shape[0]):
        m.command.start(x=b[i, 0], px=b[i, 1], y=b[i, 2], py=b[i, 3], t=0.0, pt=b[i, 4])

    m.command.observe(place='M1')

    m.command.run(turns=1)
    df = m.table.trackone.dframe()
    mr = df.loc['m1'][['x', 'px', 'y', 'py', 'pt']].values
    mr

    assert np.all(np.isclose(mr[0:4], r[:, 0:4])), True
