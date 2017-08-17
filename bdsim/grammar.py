"""A dictionnary for the translation of the MAD-X grammar."""

bdsim_syntax = {  # Do not forget the trailing ';' for each command!
    'sample': "sample, all;",
    'use': "use, period={line};",
    'beam': "beam, energy={energy}*MeV, particle=\"{particle}\";",
    'options': """option, beampipeRadius = {beampiperadius} * mm,
        apertureType ="{aperturetype}",
        beampipeThickness = {beampipethickness} * cm,
        beampipeMaterial = "{beampipematerial}",
        physicsList="em";"""
}
