"""A dictionnary for the translation of the MAD-X grammar."""

bdsim_syntax = {  # Do not forget the trailing ';' for each command!
    'sample': "sample, all;",
    'use': "use, period={line};",
    'beam': "beam, energy={energy}*MeV, particle=\"{particle}\";",
    'options': """option, beampipeRadius = {beampiperadius} * mm,
        apertureType ="{aperturetype}",
        beampipeThickness = {beampipethickness} * cm,
        beampipeMaterial = "{beampipematerial}",
        physicsList="em";""",
    'placement': """
    {line}Place: placement, sequence="{line}",
                 referenceElement = "{reference_element}",
                 referenceElementNumber = {reference_element_number},
                 x = -10*cm,
                 z = -1*m,
                 axisAngle = 1,
                 axisY = 1.,
                 angle = 0.52;
    """,
}
