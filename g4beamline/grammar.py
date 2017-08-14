"""A dictionnary for the translation of the G4BeamLine grammar."""


g4beamline_syntax = {  # Do not forget the trailing ';' for each command!

    'define_Brho' : "param Brho=2.3114",
    'define_detector' : "sample detector radius=31.5 format=ascii",
    'add_detector': "detector z={} rename={}';",

    'beam_start': "beam gaussian particle={{{{PARTICLE.lower()}}}}"
                  " beamX={} beamY={} "
                  " meanXp={} meanYp={}"
                  " meanMomentum={{{{PC*1000}}}} ",

}

