import cpymad.madx
import georges
import georges.madx
import georges_core.sequences


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


def test_madx_sequence():
    get_madx_twiss()
    madx_line = georges_core.sequences.TwissSequence(path='.',
                                                     filename='twiss.tfs',
                                                     lines=51,
                                                     with_units=True,
                                                     with_beam=True,
                                                     nparticles=100
                                                     )
    # assert georges.madx.MadX(sequence=madx_line)
