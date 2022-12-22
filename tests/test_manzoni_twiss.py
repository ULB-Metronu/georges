import os

import cpymad.madx
import numpy as np

import georges
from georges import ureg as _ureg
from georges.manzoni import Input, observers
from georges.manzoni.beam import MadXBeam

# TODO Update tests when dispersion is fixed in georges-core


def get_madx_twiss():
    m = cpymad.madx.Madx(stdout=False)
    m.input(
        """
        BEAM, PARTICLE=PROTON, ENERGY = 0.250+0.938, PARTICLE = PROTON, EX=1e-6, EY=1e-6;

        RHO:=1.35;
        KQ := +0.9;
        LCELL:=4.;
        LQ:= 0.3;
        LB:= 2.0;
        L2:=0.5*(LCELL-LB-LQ);
        L3:= 1.5;
        EANG:=10.*TWOPI/360;
        ANG := TWOPI/4;
        KQ1:= -0.9;

        D1: DRIFT, L=L3;
        D2: DRIFT, L=0.2;
        D3: DRIFT, L=0.2;
        D4: DRIFT, L=0.1;

        BD : SBEND,L=LB, ANGLE=ANG;
        MQF1 : QUADRUPOLE,L=0.5*LQ, K1=+KQ1;
        MQF2: MQF1;
        MQD : QUADRUPOLE,L=LQ, K1=-KQ1;

        !ACHROM: LINE=(MQF1, D1, MQD, D2,BD,D3, MQF2, D4);
        ACHROM: LINE=(MQF1, D1, MQD, D2, MQF2, D4);
        RING: LINE=(ACHROM);

        USE, sequence=RING;
        """,
    )
    twiss_madx = m.twiss(sequence="RING", file="twiss.tfs")
    return twiss_madx


def get_sequence():
    madx_line = georges.TwissSequence(
        path=".",
        filename="twiss.tfs",
        lines=51,
        with_units=True,
        with_beam=True,
        nparticles=500000,
        refer="exit",
    )
    return madx_line


def test_from_11_particles():
    twiss_madx = get_madx_twiss()
    madx_line = get_sequence()
    mi = Input.from_sequence(sequence=madx_line)
    mi.freeze()
    tw_observer = mi.twiss(kinematics=madx_line.kinematics)

    np.testing.assert_allclose(tw_observer["BETA11"], twiss_madx["BETX"], rtol=2e-2)
    np.testing.assert_allclose(tw_observer["BETA22"], twiss_madx["BETY"], rtol=2e-2)
    np.testing.assert_allclose(tw_observer["ALPHA11"], twiss_madx["ALFX"], rtol=2e-2)
    np.testing.assert_allclose(tw_observer["ALPHA22"], twiss_madx["ALFY"], rtol=2e-2)
    np.testing.assert_allclose(tw_observer["DISP1"] * madx_line.metadata.kinematics.beta, twiss_madx["DX"], atol=2e-2)
    np.testing.assert_allclose(tw_observer["DISP3"], twiss_madx["DY"], atol=2e-2)
    np.testing.assert_allclose(tw_observer["DISP2"], twiss_madx["DDX"] * madx_line.metadata.kinematics.beta, atol=2e-2)
    np.testing.assert_allclose(tw_observer["DISP4"], twiss_madx["DDY"], atol=2e-2)
    os.remove("twiss.tfs")


def test_from_distribution():
    twiss_madx = get_madx_twiss()
    madx_line = get_sequence()
    beam = MadXBeam(
        kinematics=madx_line.metadata.kinematics,
        distribution=georges.Distribution.from_twiss_parameters(
            n=int(1e7),
            betax=twiss_madx["BETX"][0] * _ureg.m,
            betay=twiss_madx["BETY"][0] * _ureg.m,
            alphax=twiss_madx["ALFX"][0],
            alphay=twiss_madx["ALFY"][0],
            dispx=twiss_madx["DX"][0] * _ureg.m,
            dispy=twiss_madx["DY"][0] * _ureg.m,
            dispxp=twiss_madx["DDX"][0],
            dispyp=twiss_madx["DDY"][0],
            dpprms=1e-2,
        )
        .distribution[["X", "PX", "Y", "PY", "DPP"]]
        .values,
    )
    mi = Input.from_sequence(sequence=madx_line)
    mi.freeze()
    tw_observer = mi.track(beam=beam, observers=observers.TwissObserver())
    tw_df = tw_observer.to_df()

    np.testing.assert_allclose(tw_df["BETA_OUT_X"], twiss_madx["BETX"], rtol=2e-2)
    np.testing.assert_allclose(tw_df["BETA_OUT_Y"], twiss_madx["BETY"], rtol=2e-2)
    np.testing.assert_allclose(tw_df["ALPHA_OUT_X"], twiss_madx["ALFX"], rtol=2e-2)
    np.testing.assert_allclose(tw_df["ALPHA_OUT_Y"], twiss_madx["ALFY"], rtol=2e-2)
    np.testing.assert_allclose(tw_df["DISP_OUT_X"] * madx_line.metadata.kinematics.beta, twiss_madx["DX"], atol=2e-2)
    np.testing.assert_allclose(tw_df["DISP_OUT_Y"], twiss_madx["DY"], atol=2e-2)
    np.testing.assert_allclose(tw_df["DISP_OUT_XP"], twiss_madx["DDX"] * madx_line.metadata.kinematics.beta, atol=2e-2)
    np.testing.assert_allclose(tw_df["DISP_OUT_YP"], twiss_madx["DDY"], atol=2e-2)
    os.remove("twiss.tfs")
