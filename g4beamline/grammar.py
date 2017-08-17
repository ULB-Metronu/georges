"""A dictionnary for the translation of the G4BeamLine grammar."""


g4beamline_syntax = {

    'define_brho': "param Brho={{{{BRHO}}}}",
    'define_detector': "sample detector radius=31.5 format=ascii",
    'add_detector': "detector z={} rename={}';",
    'define_physics': "physics QGSP_BIC",
    'define_world': "param worldMaterial=Vacuum",
    'keep_protons': "trackcuts keep=proton",
    'beam_start': "beam gaussian particle={{{{PARTICLE.lower()}}}}"
                  " beamX={} meanXp={} "
                  " beamY={} meanYp={} "
                  " beamZ=0.0"
                  " meanMomentum={{{{PC*1000}}}} "
                  " nEvents=1",
    'quadrupole': "genericquad {} " 
                  "fieldLength={} " 
                  "ironLength={} " 
                  "ironRadius=238 " 
                  "apertureRadius={} " 
                  "gradient={}*$Brho " 
                  "ironMaterial=Fe " 
                  "fieldMaterial=Vacuum " 
                  "ironColor=1,0,0 " 
                  "kill=1 " 
                  "openAperture=1 ", # Do not forget the whitespace
    'ccoll': "tubs {} "
            "innerRadius={} "
            "outerRadius=50 "
            "length={} "
            "material=W "
            "color=1,1,0 "
            "kill=1",
    'rcoll': "box {} " 
            "width=60 " 
            "height=80 "
            "length={} "
            "material=Ni "
            "color=1,1,0 "
            "kill=1",
}