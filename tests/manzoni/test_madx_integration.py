import pytest
import numpy as np
import cpymad.madx
from georges import ureg as _
from georges import manzoni
from georges_core import Kinematics
from georges.manzoni.integrators import MadXIntegrator


def get_madx_tracking_data_drift(x: float,
                                 px: float,
                                 y: float,
                                 py: float,
                                 pt: float,
                                 beta: float,
                                 length: float,
                                 ) -> np.ndarray:
    m = cpymad.madx.Madx(stdout=False)
    m.input(f"""
        BEAM, PARTICLE=PROTON, BETA={beta};

        SEQ: SEQUENCE, REFER=ENTRY, L={length};
          D1: DRIFT, AT=0.0, L={length};
          M1: MARKER, AT={length};
        ENDSEQUENCE;

        USE, SEQUENCE=SEQ;
        """)
    m.command.track(DELTAP=0.0, ONEPASS=True, ONETABLE=True, QUANTUM=False)
    m.command.start(x=x, px=px, y=y, py=py, t=0.0, pt=pt)
    m.command.observe(place='M1')
    m.command.run(turns=1)
    return m.table.trackone.dframe().loc['m1'][['x', 'px', 'y', 'py', 'pt']].values[0:4]


def get_madx_tracking_data_quadrupole(x: float,
                                      px: float,
                                      y: float,
                                      py: float,
                                      pt: float,
                                      beta: float,
                                      length: float,
                                      k1: float,
                                      k1s: float,
                                      tilt: float,
                                      ) -> np.ndarray:
    m = cpymad.madx.Madx(stdout=False)
    m.input(f"""
        BEAM, PARTICLE=PROTON, BETA={beta};

        SEQ: SEQUENCE, REFER=ENTRY, L={length};
          Q1: QUADRUPOLE, AT=0.0, L={length}, K1={k1}, K1S={k1s}, TILT={tilt};
          M1: MARKER, AT={length};
        ENDSEQUENCE;

        USE, SEQUENCE=SEQ;
        """)
    m.command.track(DELTAP=0.0, ONEPASS=True, ONETABLE=True, QUANTUM=False)
    m.command.start(x=x, px=px, y=y, py=py, t=0.0, pt=pt)
    m.command.observe(place='M1')
    m.command.run(turns=1)
    m.command.endtrack()
    return m.table.trackone.dframe().loc['m1'][['x', 'px', 'y', 'py', 'pt']].values[0:4]


def get_madx_tracking_data_bend(x: float,
                                px: float,
                                y: float,
                                py: float,
                                pt: float,
                                beta: float,
                                length: float,
                                angle: float,
                                k1: float,
                                e1: float = 0.0,
                                e2: float = 0.0,
                                hgap: float = 0.0,
                                fint: float = 0.0,
                                fintx: float = 0.0,
                                element: str = 'SBEND',
                                ) -> np.ndarray:
    if element == 'RBEND' and angle > 1e-8:
        arc_length = length * angle / (2.0 * np.sin(angle / 2.0))
    else:
        arc_length = length
    m = cpymad.madx.Madx(stdout=False)
    m.input(f"""
    BEAM, PARTICLE=PROTON, BETA={beta};

    SEQ: SEQUENCE, REFER=ENTRY, L={arc_length};
      B1: {element}, AT=0.0, THICK=true, L={length}, ANGLE={angle}, K1={k1}, E1={e1}, E2={e2}, HGAP={hgap}, FINT={fint}, FINTX={fintx};
      M1: MARKER, AT={arc_length};
    ENDSEQUENCE;

    USE, SEQUENCE=SEQ;
        """)
    m.command.track(DELTAP=0.0, ONEPASS=True, ONETABLE=True, QUANTUM=False)
    m.command.start(x=x, px=px, y=y, py=py, t=0.0, pt=pt)
    m.command.observe(place='M1')
    m.command.run(turns=1)
    m.command.endtrack()
    return m.table.trackone.dframe().loc['m1'][['x', 'px', 'y', 'py', 'pt']].values[0:4]


def get_madx_tracking_data_dipedge(x: float,
                                   px: float,
                                   y: float,
                                   py: float,
                                   pt: float,
                                   beta: float,
                                   h: float,
                                   e1: float,
                                   hgap: float,
                                   fint: float,
                                   ) -> np.ndarray:
    m = cpymad.madx.Madx(stdout=False)
    m.input(f"""
    BEAM, PARTICLE=PROTON, BETA={beta};

    SEQ: SEQUENCE, REFER=ENTRY, L=1.0;
      D1: DRIFT, AT=0.0, L=0.5;
      DE1: DIPEDGE, AT=0.5, H={h}, E1={e1}, HGAP={hgap}, FINT={fint};
      D2: DRIFT, AT=0.5, L=0.5;
      M1: MARKER, AT={1.0};
    ENDSEQUENCE;

    USE, SEQUENCE=SEQ;
        """)
    m.command.track(DELTAP=0.0, ONEPASS=True, ONETABLE=True, QUANTUM=False)
    m.command.start(x=x, px=px, y=y, py=py, t=0.0, pt=pt)
    m.command.observe(place='M1')
    m.command.run(turns=1)
    m.command.endtrack()
    return m.table.trackone.dframe().loc['m1'][['x', 'px', 'y', 'py', 'pt']].values[0:4]


def get_madx_tracking_data_kicker(x: float,
                                  px: float,
                                  y: float,
                                  py: float,
                                  pt: float,
                                  beta: float,
                                  length: float,
                                  tilt: float,
                                  hkick: float,
                                  vkick: float) -> np.ndarray:
    m = cpymad.madx.Madx(stdout=False)
    m.input(f"""
    BEAM, PARTICLE=PROTON, BETA={beta};

    SEQ: SEQUENCE, REFER=ENTRY, L=1.0;
      K1: KICKER, AT=0.0, L={length}, TILT={tilt}, HKICK={hkick}, VKICK={vkick};
      M1: MARKER, AT={length};
    ENDSEQUENCE;

    USE, SEQUENCE=SEQ;
        """)
    m.command.track(DELTAP=0.0, ONEPASS=True, ONETABLE=True, QUANTUM=False)
    m.command.start(x=x, px=px, y=y, py=py, t=0.0, pt=pt)
    m.command.observe(place='M1')
    m.command.run(turns=1)
    m.command.endtrack()
    return m.table.trackone.dframe().loc['m1'][['x', 'px', 'y', 'py', 'pt']].values[0:4]


@pytest.mark.parametrize("x, px, y, py, pt, beta, length", [
    (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0),
    (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 10.0),
    (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 100.0),
    (0.1, 0.0, 0.0, 0.0, 0.0, 0.9, 1.0),
    (0.1, 0.0, 0.0, 0.0, 0.0, 0.9, 10.0),
    (0.1, 0.0, 0.0, 0.0, 0.0, 0.9, 100.0),
    (1.0, 0.0, 0.0, 0.0, 0.0, 0.5, 100.0),
    (0.0, 0.1, 0.0, 0.0, 0.0, 0.5, 100.0),
    (0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 100.0),
    (0.0, 0.0, 0.0, 0.1, 0.0, 0.5, 100.0),
    (0.0, 0.0, 0.0, 0.0, 1.0, 0.5, 100.0),
    (1.0, 0.1, 0.0, 0.0, 0.0, 0.5, 100.0),
    (1.0, 0.1, 1.0, 0.0, 0.0, 0.5, 100.0),
    (1.0, 0.1, 1.0, 0.0, 0.0, 0.5, 100.0),
    (1.0, 0.1, 1.0, 0.1, 0.1, 0.5, 51.0),
    (1.0, 0.1, 1.0, 0.1, 0.1, 0.5, 11.0),
    (0.0, -0.01, 0.0, 0.0, 0.01, 0.999999999999, 1.0),
    (1.0, 0.0, 1.0, 0.0, 0.1, 0.99, 51.0),
    (1.0, 0.1, 1.0, 0.1, 0.1, 0.99, 11.0),
    (0.0, 0.1, 0.0, 0.1, 0.1, 0.99, 101.0),
    (1.0, 0.0, 1.0, 0.0, 0.1, 0.5, 100.0),
    (1.0, 0.0, 1.0, 0.0, 1.0, 0.5, 100.0),
    (1.0, 0.1, 1.0, 0.1, 0.0, 0.5, 100.0),
])
def test_madx_drift(x, px, y, py, pt, beta, length):
    k = Kinematics(beta)
    mi = manzoni.Input(sequence=[
        manzoni.Drift('B1', integrator=MadXIntegrator, L=length * _.m),
    ]).freeze()
    deltap = np.sqrt(pt ** 2 + 2 * pt / beta + 1) - 1
    distribution = np.array([
        [x, px, y, py, deltap, pt],
    ])
    beam = manzoni.Beam(kinematics=k, distribution=distribution)
    observer = manzoni.BeamObserver()
    manzoni.track(beamline=mi, beam=beam, observer=observer)
    r = observer.to_df().iloc[-1]['BEAM_OUT']

    assert np.all(np.isclose(get_madx_tracking_data_drift(
        x, px, y, py, pt, beta, length
    ), r[:, 0:4]))


@pytest.mark.parametrize(
    "x, px, y, py, pt, beta, length, k1, k1s, tilt", [
        # Field free quadrupole
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 0.0, 0.0, 0.0),

        # Regular quadrupole
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, -1.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 10.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, -10.0, 0.0, 0.0),

        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, -1.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 10.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, -10.0, 0.0, 0.0),

        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, -1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, 10.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, -10.0, 0.0, 0.0),

        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, -1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 10.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, -10.0, 0.0, 0.0),

        (0.1, 0.0, 0.0, 0.0, 0.2, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.2, 0.5, 1.0, 1.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.2, 0.5, 1.0, -1.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.2, 0.5, 1.0, 10.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.2, 0.5, 1.0, -10.0, 0.0, 0.0),

        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, -1.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 10.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, -10.0, 0.0, 0.0),

        (0.1, 0.1, 0.1, 0.0, 0.2, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.2, 0.5, 1.0, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.2, 0.5, 1.0, -1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.2, 0.5, 1.0, 10.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.2, 0.5, 1.0, -10.0, 0.0, 0.0),

        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, -1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 10.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, -10.0, 0.0, 0.0),

        # Quadrupole with pure skew component
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 1.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, -1.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 10.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, -10.0, 0.0),

        # Quadrupole with mixed normal and skew components
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, -1.0, -1.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 10.0, 10.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, -10.0, -10.0, 0.0),

        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, -1.0, -1.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 10.0, 10.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, -10.0, -10.0, 0.0),

        # Quadrupole with tilt
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.1),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 0.0, 0.1),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, -1.0, 0.0, 0.1),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 10.0, 0.0, 0.1),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, -10.0, 0.0, 0.1),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, -0.1),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 0.0, -0.1),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, -1.0, 0.0, -0.1),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 10.0, 0.0, -0.1),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, -10.0, 0.0, -0.1),
    ])
def test_madx_quadrupole(x, px, y, py, pt, beta, length, k1, k1s, tilt):
    k = Kinematics(beta)
    mi = manzoni.Input(sequence=[
        manzoni.Quadrupole('B1',
                           integrator=MadXIntegrator,
                           L=length * _.m,
                           K1=k1 * _.m**-2,
                           K1S=k1s * _.m**-2,
                           TILT=tilt * _.radian
                           ),
    ]).freeze()
    deltap = np.sqrt(pt ** 2 + 2 * pt / beta + 1) - 1
    distribution = np.array([
        [x, px, y, py, deltap, pt],
    ])
    beam = manzoni.Beam(kinematics=k, distribution=distribution)
    observer = manzoni.BeamObserver()
    manzoni.track(beamline=mi, beam=beam, observer=observer)
    r = observer.to_df().iloc[-1]['BEAM_OUT']

    assert np.all(np.isclose(get_madx_tracking_data_quadrupole(
        x, px, y, py, pt, beta, length, k1, k1s, tilt
    ), r[:, 0:4]))


@pytest.mark.parametrize(
    "x, px, y, py, pt, beta, length, angle, k1", [
        # Regular bend
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.1, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 0.0, 0.0),

        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.1, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 0.1, 0.0),

        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.0),
        (0.0, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.1, 0.0, 0.1, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 0.1, 1.0),
    ])
def test_madx_sbend(x, px, y, py, pt, beta, length, angle, k1):
    k = Kinematics(beta)
    mi = manzoni.Input(sequence=[
        manzoni.SBend('B1',
                      integrator=MadXIntegrator,
                      L=length * _.m,
                      K1=k1 * _.m ** -2,
                      ANGLE=angle * _.radian,
                      ),
    ]).freeze()
    deltap = np.sqrt(pt ** 2 + 2 * pt / beta + 1) - 1
    distribution = np.array([
        [x, px, y, py, deltap, pt],
    ])
    beam = manzoni.Beam(kinematics=k, distribution=distribution)
    observer = manzoni.BeamObserver()
    manzoni.track(beamline=mi, beam=beam, observer=observer)
    r = observer.to_df().iloc[-1]['BEAM_OUT']

    assert np.all(np.isclose(get_madx_tracking_data_bend(
        x, px, y, py, pt, beta, length, angle, k1
    ), r[:, 0:4]))


@pytest.mark.parametrize(
    "x, px, y, py, pt, beta, length, angle, k1", [
        # Regular bend
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.1, 0.5, 1.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 0.0, 0.0),

        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.0, 0.1, 0.5, 1.0, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 0.1, 0.0),

        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.0),
        (0.0, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.1, 0.0, 0.1, 0.5, 1.0, 0.1, 1.0),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 0.1, 1.0),
    ])
def test_madx_rbend(x, px, y, py, pt, beta, length, angle, k1):
    k = Kinematics(beta)
    mi = manzoni.Input(sequence=[
        manzoni.RBend('B1',
                      integrator=MadXIntegrator,
                      L=length * _.m,
                      K1=k1 * _.m ** -2,
                      ANGLE=angle * _.radian,
                      ),
    ]).freeze()
    deltap = np.sqrt(pt ** 2 + 2 * pt / beta + 1) - 1
    distribution = np.array([
        [x, px, y, py, deltap, pt],
    ])
    beam = manzoni.Beam(kinematics=k, distribution=distribution)
    observer = manzoni.BeamObserver()
    manzoni.track(beamline=mi, beam=beam, observer=observer)
    r = observer.to_df().iloc[-1]['BEAM_OUT']

    assert np.all(np.isclose(get_madx_tracking_data_bend(
        x, px, y, py, pt, beta, length, angle, k1, element='RBEND'
    ), r[:, 0:4]))


@pytest.mark.parametrize(
    "x, px, y, py, pt, beta, length, angle, k1, e1, e2, hgap, fint, fintx", [
        # X    PX   Y    PY   PT  BETA  L   ANG   K1   E1  E2   HGAP FINT FINTX
        # No bending angle
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        # Only E1
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        # Only E2
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.1, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.1, 0.0, 0.1),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.1, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.1, 0.0, 0.1),
        # E1 and E2
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.1, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.1, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, -0.1, -0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, -0.1, -0.1, 0.0, 0.0, 0.0),
        # With field integral (entrance) only
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        # With field integral (exit) only
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.0, 0.5),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.0, 0.5),
    ])
def test_madx_sbend_fringe(x, px, y, py, pt, beta, length, angle, k1, e1, e2, hgap, fint, fintx):
    k = Kinematics(beta)
    mi = manzoni.Input(sequence=[
        manzoni.SBend('B1',
                      integrator=MadXIntegrator,
                      L=length * _.m,
                      K1=k1 * _.m ** -2,
                      ANGLE=angle * _.radian,
                      E1=e1 * _.radian,
                      E2=e2 * _.radian,
                      HGAP=hgap * _.m,
                      FINT=fint,
                      FINTX=fintx,
                      ),
    ]).freeze()
    deltap = np.sqrt(pt ** 2 + 2 * pt / beta + 1) - 1
    distribution = np.array([
        [x, px, y, py, deltap, pt],
    ])
    beam = manzoni.Beam(kinematics=k, distribution=distribution)
    observer = manzoni.BeamObserver()
    manzoni.track(beamline=mi, beam=beam, observer=observer)
    r = observer.to_df().iloc[-1]['BEAM_OUT']

    assert np.all(np.isclose(get_madx_tracking_data_bend(
        x, px, y, py, pt, beta, length, angle, k1, e1, e2, hgap, fint, fintx
    ), r[:, 0:4]))


@pytest.mark.parametrize(
    "x, px, y, py, pt, beta, length, angle, k1, e1, e2, hgap, fint, fintx", [
        # X    PX   Y    PY   PT  BETA  L   ANG   K1   E1  E2   HGAP FINT FINTX
        # No bending angle
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        # Only E1
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        # Only E2
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.1, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.1, 0.0, 0.1),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.1, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.1, 0.0, 0.1),
        # E1 and E2
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.1, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.1, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, -0.1, -0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, -0.1, -0.1, 0.0, 0.0, 0.0),
        # With field integral (entrance) only
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        # With field integral (exit) only
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.0, 0.5),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.0, 0.5),
    ])
def test_madx_sbend_fringe(x, px, y, py, pt, beta, length, angle, k1, e1, e2, hgap, fint, fintx):
    k = Kinematics(beta)
    mi_sbend = manzoni.Input(sequence=[
        manzoni.SBend('B1',
                      integrator=MadXIntegrator,
                      L=length * _.m,
                      K1=k1 * _.m ** -2,
                      ANGLE=angle * _.radian,
                      E1=e1 * _.radian,
                      E2=e2 * _.radian,
                      HGAP=hgap * _.m,
                      FINT=fint,
                      FINTX=fintx,
                      ),
    ]).freeze()
    deltap = np.sqrt(pt ** 2 + 2 * pt / beta + 1) - 1
    distribution = np.array([
        [x, px, y, py, deltap, pt],
    ])
    beam = manzoni.Beam(kinematics=k, distribution=distribution)
    observer = manzoni.BeamObserver()
    manzoni.track(beamline=mi_sbend, beam=beam, observer=observer)
    r_sbend = observer.to_df().iloc[-1]['BEAM_OUT']

    assert np.all(np.isclose(get_madx_tracking_data_bend(
        x, px, y, py, pt, beta, length, angle, k1, e1, e2, hgap, fint, fintx, 'SBEND'
    ), r_sbend[:, 0:4]))


@pytest.mark.parametrize(
    "x, px, y, py, pt, beta, length, angle, k1, e1, e2, hgap, fint, fintx", [
        # X  PX   Y     PY   PT  BETA  L   ANG   K1   E1  E2   HGAP FINT FINTX
        # No bending angle
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.0, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        # Only E1
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.1, 1.3, 0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.1, 1.3, -0.1, 0.0, 0.0, 0.0, 0.0),
        # Only E2
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.1, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.1, 0.1, 0.0, 0.1),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.1, 0.1, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, -0.1, 0.1, 0.0, 0.1),
        # E1 and E2
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.1, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.1, 0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, -0.1, -0.1, 0.0, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, -0.1, -0.1, 0.0, 0.0, 0.0),
        # With field integral (entrance) only
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.5, 0.0),
        # With field integral (exit) only
        (0.1, 0.1, 0.1, 0.1, 0.0, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.0, 0.5),
        (0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 1.0, 0.1, 1.3, 0.0, 0.0, 0.1, 0.0, 0.5),
    ])
def test_madx_rbend_fringe(x, px, y, py, pt, beta, length, angle, k1, e1, e2, hgap, fint, fintx):
    k = Kinematics(beta)
    mi_rbend = manzoni.Input(sequence=[
        manzoni.RBend('B1',
                      integrator=MadXIntegrator,
                      L=length * _.m,
                      K1=k1 * _.m ** -2,
                      ANGLE=angle * _.radian,
                      E1=e1 * _.radian,
                      E2=e2 * _.radian,
                      HGAP=hgap * _.m,
                      FINT=fint,
                      FINTX=fintx,
                      ),
    ]).freeze()
    deltap = np.sqrt(pt ** 2 + 2 * pt / beta + 1) - 1
    distribution = np.array([
        [x, px, y, py, deltap, pt],
    ])
    beam = manzoni.Beam(kinematics=k, distribution=distribution)
    observer = manzoni.BeamObserver()
    manzoni.track(beamline=mi_rbend, beam=beam, observer=observer)
    r_rbend = observer.to_df().iloc[-1]['BEAM_OUT']

    assert np.all(np.isclose(get_madx_tracking_data_bend(
        x, px, y, py, pt, beta, length, angle, k1, e1, e2, hgap, fint, fintx, 'RBEND'
    ), r_rbend[:, 0:4]))


@pytest.mark.parametrize(
    "x, px, y, py, pt, beta, h, e1, hgap, fint", [
        # X   PX    Y    PY   PT  BETA H    E1  HGAP  FINT
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, 0.1, 0.1, 0.5),
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, 0.0, 0.1, 0.5),
        (0.0, 0.0, 0.1, 0.1, 0.1, 0.5, 1.0, 0.1, 0.0, 0.0),
        (0.0, 0.0, 0.1, 0.1, 0.1, 0.5, 1.0, 0.1, 0.1, 0.5),
        (0.0, 0.0, 0.1, 0.1, 0.1, 0.5, 1.0, 0.0, 0.1, 0.5),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 0.1, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 0.1, 0.1, 0.5),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 0.0, 0.1, 0.5),
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, -0.1, 0.0, 0.0),
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, -0.1, 0.1, 0.5),
        (0.1, 0.1, 0.0, 0.0, 0.1, 0.5, 1.0, -0.0, 0.1, 0.5),
        (0.0, 0.0, 0.1, 0.1, 0.1, 0.5, 1.0, -0.1, 0.0, 0.0),
        (0.0, 0.0, 0.1, 0.1, 0.1, 0.5, 1.0, -0.1, 0.1, 0.5),
        (0.0, 0.0, 0.1, 0.1, 0.1, 0.5, 1.0, -0.0, 0.1, 0.5),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, -0.1, 0.0, 0.0),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, -0.1, 0.1, 0.5),
        (0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, -0.0, 0.1, 0.5),
    ])
def test_madx_dipedge(x, px, y, py, pt, beta, h, e1, hgap, fint):
    k = Kinematics(beta)
    mi = manzoni.Input(sequence=[
        manzoni.Drift('D1', L=0.5 * _.m),
        manzoni.DipEdge('DE1',
                        integrator=MadXIntegrator,
                        H=h * _.m ** -1,
                        E1=e1 * _.radian,
                        HGAP=hgap * _.m,
                        FINT=fint,
                        ),
        manzoni.Drift('D2', L=0.5 * _.m),
    ]).freeze()
    deltap = np.sqrt(pt ** 2 + 2 * pt / beta + 1) - 1
    distribution = np.array([
        [x, px, y, py, deltap, pt],
    ])
    beam = manzoni.Beam(kinematics=k, distribution=distribution)
    observer = manzoni.BeamObserver()
    manzoni.track(beamline=mi, beam=beam, observer=observer)
    r = observer.to_df().iloc[-1]['BEAM_OUT']

    assert np.all(np.isclose(get_madx_tracking_data_dipedge(
        x, px, y, py, pt, beta, h, e1, hgap, fint
    ), r[:, 0:4]))


@pytest.mark.parametrize(
    "x, px, y, py, pt, beta, length, tilt, hkick, vkick", [
        # X   PX    Y    PY   PT  BETA L  TILT HKICK VKICK
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (1.0, 0.0, 1.0, 0.0, 1.0, 0.5, 0.0, 0.0, 0.0, 0.0),
        (1.0, 0.0, 1.0, 0.0, 1.0, 0.5, 1.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.5, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.5),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.5, 0.5),
        (1.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0),
        (1.0, 0.0, 1.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.5, 0.0),
        (1.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5),
        (1.0, 0.0, 1.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.5),
        (1.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5),
        (1.0, 0.0, 1.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.5, 0.5),
        # With tilt
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.0, 0.0),
        (1.0, 0.1, 1.0, 0.0, 1.0, 0.5, 0.0, 0.1, 0.0, 0.0),
        (1.0, 0.1, 1.0, 0.0, 1.0, 0.5, 1.0, 0.1, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.1, 0.5, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.5, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.0, 0.5),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.1, 0.5, 0.5),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.5, 0.5),
        (1.0, 0.1, 1.0, 0.0, 0.0, 0.5, 0.0, 0.1, 0.5, 0.0),
        (1.0, 0.1, 1.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.5, 0.0),
        (1.0, 0.1, 1.0, 0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5),
        (1.0, 0.1, 1.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.0, 0.5),
        (1.0, 0.1, 1.0, 0.0, 0.0, 0.5, 0.0, 0.1, 0.5, 0.5),
        (1.0, 0.1, 1.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.5, 0.5),
    ])
def test_madx_kicker(x, px, y, py, pt, beta, length, tilt, hkick, vkick):
    k = Kinematics(beta)
    mi = manzoni.Input(sequence=[
        manzoni.Kicker('K1',
                       integrator=MadXIntegrator,
                       L=length * _.m,
                       TILT=tilt * _.radian,
                       HKICK=hkick * _.radian,
                       VKICK=vkick * _.radian,
                       ),
    ]).freeze()
    deltap = np.sqrt(pt ** 2 + 2 * pt / beta + 1) - 1
    distribution = np.array([
        [x, px, y, py, deltap, pt],
    ])
    beam = manzoni.Beam(kinematics=k, distribution=distribution)
    observer = manzoni.BeamObserver()
    manzoni.track(beamline=mi, beam=beam, observer=observer)
    r = observer.to_df().iloc[-1]['BEAM_OUT']

    assert np.all(np.isclose(get_madx_tracking_data_kicker(
        x, px, y, py, pt, beta, length, tilt, hkick, vkick
    ), r[:, 0:4]))
