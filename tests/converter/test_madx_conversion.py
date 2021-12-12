import os
import cpymad.madx
import numpy as np
import pandas as pd
import georges
from georges.manzoni import Input, Beam
from georges.manzoni import observers


def get_madx_twiss():
    m = cpymad.madx.Madx(stdout=False)
    m.input(f"""
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
    
        OO : DRIFT,L=L2;
        BD : SBEND,L=LB, ANGLE=ANG, E1=EANG,E2=EANG, K2=0.;
        MQ : QUADRUPOLE,L=LQ, K1=KQ;
        MQF : QUADRUPOLE,L=LQ, K1=+KQ1;
        MQD : QUADRUPOLE,L=LQ, K1=-KQ1;
        OBS1: MARKER;
        OBS2: MARKER;
    
        ACHROM: LINE=(MQD,D2,MQF,D2,BD,OO,MQ,OO,BD,D2,MQF,D2,MQD,OBS1);
        RING: LINE=(ACHROM, D1,ACHROM,D1);
    
        USE, sequence=RING;
        """)
    m.command.twiss(sequence='RING', FILE='twiss.tfs')


def test_madx_conversion():
    get_madx_twiss()
    madx_line = georges.TwissSequence(path='.',
                                      filename='twiss.tfs',
                                      with_units=True,
                                      with_beam=True,
                                      nparticles=500000,
                                      refer='exit'
                                      )
    beam = Beam(kinematics=madx_line.metadata.kinematics,
                distribution=madx_line.metadata.data.values
                )

    mi = Input.from_sequence(sequence=madx_line)
    tw_observer = mi.track(beam=beam, observers=observers.TwissObserver())

    tw_df = tw_observer.to_df()
    md_df = madx_line.df
    idx = tw_df.index.values
    for i, k in enumerate(idx):
        idx[i] += f"{i}"
    tw_df.set_index(idx, inplace=True)
    md_df.set_index(tw_df.index, inplace=True)

    df = pd.merge(tw_df, md_df, left_index=True, right_index=True)[['BETX', 'BETY', 'ALFX', 'ALFY',
                                                                    'BETA_OUT_X', 'BETA_OUT_Y', 'ALPHA_OUT_X',
                                                                    'ALPHA_OUT_Y']]
    np.testing.assert_allclose(df['BETX'], df['BETA_OUT_X'], rtol=2e-2)
    np.testing.assert_allclose(df['BETY'], df['BETA_OUT_Y'], rtol=2e-2)
    np.testing.assert_allclose(df['ALFX'], df['ALPHA_OUT_X'], rtol=1)
    np.testing.assert_allclose(df['ALFY'], df['ALPHA_OUT_Y'], rtol=1)
    os.remove("twiss.tfs")