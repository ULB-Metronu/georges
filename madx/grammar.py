"""A dictionnary for the translation of the MAD-X grammar."""

madx_syntax = {  # Do not forget the trailing ';' for each command!
    'beam': "BEAM, PARTICLE={{{{PARTICLE}}}}, "
            "PC={{{{PC/1000.0}}}}, "
            "EX={{{{ EMITX or '1e-9' }}}}, "
            "EY={{{{ EMITY or '1e-9' }}}};",
    'show_beam': "SHOW, BEAM;",
    'call_file': "CALL, FILE='{}';",
    'use_sequence': "USE, SEQUENCE={};",
    'save_beta': "SAVEBETA, LABEL={}, PLACE={};",
    'misalign_option':"EOPTION, ADD=False;",
    'mad_misalign_setup': "SELECT, FLAG=ERROR,CLASS={};\n"
                          "EALIGN, "
                          "DX={},"
                          "DY={},"
                          "DS={},"
                          "DPHI={},"
                          "DTHETA={};",
    'makethin': "MAKETHIN, sequence={}, style={};",
    'twiss_beamline': "TWISS, "
                      "BETX={{{{ BETAX or '1.0' }}}},"
                      "ALFX={{{{ ALPHAX or '0.0' }}}},"
                      "MUX=0.0,"
                      "BETY={{{{ BETAY or '1.0' }}}},"
                      "ALFY={{{{ ALPHAY or '0.0' }}}},"
                      "MUY=0.0,"
                      "DX=0.0,"
                      "DPX=0.0,"
                      "DY=0.0,"
                      "DPY=0.0,"
                      "X=0.0,"
                      "PX=0.0,"
                      "Y=0.0,"
                      "PY=0.0,"
                      "T=0.0,"
                      "PT=0.0,"
                      "DELTAP={{{{ DELTAP or '0.0' }}}},"
                      "FILE={},"
                      "SECTORMAP=true" # NO coma
                      "{};",  # Note the optional args
    'track_beamline': "TRACK,"
                      "DELTAP={{{{ DELTAP or '0.0' }}}},"
                      "ONEPASS=true,"
                      "APERTURE=true,"
                      "DUMP=true,"
                      "ONETABLE=true,"
                      "FILE=tracking.outx;",
    'ptc_twiss': "PTC_TWISS,ICASE=56,"
                 "DELTAP={{{{ DELTAP or '0.0' }}}},"
                 "FILE={},"
                 "CLOSED_ORBIT=true,"
                 "NO=4,"
                 "DELTAP_DEPENDENCY=true,"
                 "SLICE_MAGNETS=true;",
    'ptc_twiss_beamline': "PTC_TWISS,ICASE=56,"
                          "NO=2,"
                          "DELTAP={{{{ DELTAP or '0.0' }}}},"
                          "FILE={},"
                          "BETX={{{{ BETAX or '1.0' }}}},"
                          "ALFX={{{{ ALPHAX or '0.0' }}}},"
                          "MUX={{{{ MUX or '0.0' }}}},"
                          "BETY={{{{ BETAY or '1.0' }}}},"
                          "ALFY={{{{ ALPHAY or '0.0' }}}},"
                          "MUY={{{{ MUY or '0.0' }}}},"
                          "DX={{{{ DISP1 or '0.0' }}}},"
                          "DPX={{{{ DISP2 or '0.0' }}}},"
                          "DY={{{{ DISP3 or '0.0' }}}},"
                          "DPY={{{{ DISP4 or '0.0' }}}},"
                          "X={{{{ X or '0.0' }}}},"
                          "PX={{{{ PX or '0.0' }}}},"
                          "Y={{{{ Y or '0.0' }}}},"
                          "PY={{{{ PY or '0.0' }}}},"
                          "T=0.0,"
                          "BETZ=1.0,"
                          "ALFZ=0.0,"
                          "MUZ=0.0,"
                          "PT=0.0,"
                          "SLICE_MAGNETS=true;",
    'run_track_beamline': "RUN, TURNS=1, MAXAPER={0.1, 0.01, 0.1, 0.01, 1.0, 0.1};",  # Beamline so OK to hardcode TURNS=1
    'start_particle': "START, X={}, PX={}, Y={}, PY={}, T=0.0, PT={};",
    'observe': "OBSERVE, PLACE={};",
    'end_track': 'ENDTRACK;',
    'stop': "STOP;",
    'rbarc': "OPTION, RBARC=false;",
    'select_columns': "SELECT, FLAG={}, COLUMN={};",
    'eager_variable': "{} = {{{{ {} }}}};",  # Oops
    'lazy_variable': "{} := {{{{ {} }}}};",  # Oops
    'ptc_create_universe': "PTC_CREATE_UNIVERSE;",
    'ptc_create_layout': "PTC_CREATE_LAYOUT, TIME={}, MODEL={}, METHOD={}, NST={}, EXACT={};\nPTC_SETSWITCH, FRINGE={};",
    'ptc_misalign': "PTC_ALIGN;",
    'ptc_align': "PTC_ALIGN;",
    'ptc_end': "PTC_END;",
    'ptc_observe': "PTC_OBSERVE, PLACE={};",
    'ptc_start': "PTC_START, X={}, PX={}, Y={}, PY={}, T=0.0, PT={};",
    'ptc_track': "PTC_TRACK, ICASE={},"
                 "DELTAP={},"
                 "CLOSED_ORBIT={},"
                 "ELEMENT_BY_ELEMENT={},"
                 "TURNS={},"
                 "MAXAPER={{1.0, 100, 1.0,100,100,100}},"
                 "DUMP={},"
                 "ONETABLE={},"
                 "FILE={},"
                 "EXTENSION={};",
    'ptc_track_end': "PTC_TRACK_END;",
    'match_line': "MATCH,SEQUENCE={},"
             "BETX={},ALFX={},MUX={},"
             "BETY={},ALFY={},MUY={},"
             "X={},PX={},Y={},PY={},"
             "DX={},DY={},DPX={},DPY={},"
             "DELTAP={};",
    'match_vary': "VARY,NAME={}, LOWER={}, UPPER={};",
    'match_constraint': "CONSTRAINT,RANGE='{}',{};",
    'match_jacobian': "JACOBIAN, CALLS=5000, TOLERANCE=1E-6, REPEAT=1,"
                      "STRATEGY=1;",
    'match_lmdif': "LMDIF, CALLS=1000, TOLERANCE=1E-6;",
    'match_migrad': "MIGRAD, CALLS=1000, TOLERANCE=1E-6, STRATEGY=1;",
    'end_match': "ENDMATCH;",
    'survey': "SURVEY, file=survey.out;",
}
