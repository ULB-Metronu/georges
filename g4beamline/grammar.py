"""A dictionnary for the translation of the G4BeamLine grammar."""


g4beamline_syntax = {

    'define_detector': "sample detector radius=31.5 format=ascii",
    'define_physics': "physics QGSP_BIC",
    'define_world': "param worldMaterial=Vacuum",
    'keep_protons': "trackcuts keep=proton",
    'start_command': "start x=0 y=0 z=0 initialZ=0 radiusCut=31.5",  # Define the equivalent of MaxAper in G4Beamline
    'beam_input': "beam ascii file=input_beam.dat nEvents={}",

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
                  "gradient={} " 
                  "ironMaterial=Fe " 
                  "fieldMaterial=Vacuum " 
                  "ironColor=1,0,0 " 
                  "kill=1 " 
                  "openAperture=1 ",  # Do not forget the whitespace
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

    'tesselatedsolids': "tessellatedsolid {} kill=0 material={} file={}",
    'place_solids': "place {} x={} y={} z={} rotation=X{},Y{},Z{}",
    'virtual_det': " virtualdetector det radius=31.5 innerRadius=0 \n place det z={} rename=det_{} format=ascii \n"
}
