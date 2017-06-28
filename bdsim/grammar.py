"""A dictionnary for the translation of the MAD-X grammar."""

bdsim_syntax = {  # Do not forget the trailing ';' for each command!
    'sample': "sample, all;",
    'beam': "beam, particle={}",
    'options': """option, beampipeRadius = {} * mm,
        apertureType ="{}",
        beampipeThickness = {} * cm,
        beampipeMaterial = "{}",
        physicsList="em";"""
}
