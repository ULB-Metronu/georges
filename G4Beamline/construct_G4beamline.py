# bl.line.to_csv('test.csv')
# print(bl.line)

with open('Beamline.g4bl', 'w') as file:
    BRho = beamphys.momentum_to_brho(bl.context['PC'])
    fringe = '0'
    if (fringe):
        fringe = '1'
    for name, element in bl.line.iterrows():
        if (element['TYPE'] == "QUADRUPOLE"):
            quadrupole_cmd = ('genericquad ' + name + ' '
                                                      'fieldLength=' + str(element['LENGTH'] * 1000) + ' '
                                                                                                       'ironLength=' + str(
                element['LENGTH'] * 1000) + ' '
                                            'ironRadius=238 ' +
                              'apertureRadius=' + str(element['APERTURE'] * 1000) + ' '
                                                                                    'gradient=' + str(
                bl.context[element['CIRCUIT']] * BRho) + ' '
                                                         'ironMaterial=Fe ' +
                              'fieldMaterial=Vacuum ' +
                              'ironColor=1,0,0 ' +
                              'kill=1 ' +
                              'fringe=' + fringe + ' '
                                                   'openAperture=1'
                              )
            place_cmd = ('place ' + name + ' '
                                           'z=' + str(element['AT_CENTER'] * 1000)
                         )

            file.write(quadrupole_cmd)
            file.write('\n')
            file.write(place_cmd)
            file.write('\n')

