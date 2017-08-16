"""A dictionnary for the translation of the G4BeamLine grammar."""


g4beamline_syntax = {  # Do not forget the trailing ';' for each command!

    'define_brho': "param Brho={{{{BRHO}}}}",
    'define_detector': "sample detector radius=31.5 format=ascii",
    'add_detector': "detector z={} rename={}';",
    'define_physics': "physics QGSP_BIC",
    'define_world': "param worldMaterial=Vacuum",
    'keep_protons': "trackcuts keep=proton",
    'beam_start': "beam gaussian particle={{{{PARTICLE.lower()}}}}"
                  " beamX={} beamY={} beamZ=0.0"
                  " meanXp={} meanYp={}"
                  " meanMomentum={{{{PC*1000}}}} "
                  " nEvents=1",

}

