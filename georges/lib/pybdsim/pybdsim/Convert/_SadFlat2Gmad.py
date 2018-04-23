import numpy as _np
import optparse as _op
from collections import OrderedDict

import pysad

from .. import Builder 
from .. import Beam
from .. import Options
from .. import XSecBias

def SadFlat2GMad(inputFileName,
                 outputFileName,
                 istart=0,
                 samplers='all',
                 options = True,
                 enableSextupoles=False) :
    r = pysad.Reader.Flat(inputFileName)

    # Need nominal energy for SR calculations
    energy = 1.282
    # Emittance
    gemit    = _np.zeros(2)
    gemit[0] = 1.2e-9
    gemit[1] = 12e-12

    # Energy spread
    esprd    = 0.0006

    # particle and flip required
    particle = 'e-'
    flip     = -1.0

    # Twiss parameters
    beam = Beam.Beam(particle,energy,'gausstwiss')
    beam._SetEmittanceX(gemit[0], 'm')
    beam._SetEmittanceY(gemit[1], 'm')
    beam._SetBetaX(r.row(0)['BX'])
    beam._SetBetaY(r.row(0)['BY'])
    beam._SetAlphaX(r.row(0)['AX'])
    beam._SetAlphaY(r.row(0)['AY'])
    beam._SetSigmaE(esprd)

    # Default beam pipe radius
    defaultBpRadius = 0.05

    # create machine instance
    a = Builder.Machine(sr=False, energy0=float(energy))    

    # add the beam
    a.AddBeam(beam)


    # create options
    if options:
        o = Options.ElectronColliderOptions()
        o.SetBuildTunnel(False)
        o.SetBuildTunnelFloor(False)
        o.SetMagnetGeometryType("none")

        process = 'em'
        o.SetPhysicsList(process)
        a.AddOptions(o)

    # create name dictionary
    nameDict = {}
    unknDict = {}

    # iterate through objects and build machine
    for i in range(istart, len(r.column('NAME')), 1):

        row = r.row(i)

        # check name dictionary
        name = row['NAME']
        try:
            eCount = nameDict[name]
        except KeyError:
            nameDict[name] = 0
            eCount = nameDict[name]

        # check if name starts with a number
        prepend = ''
        if name[0].isdigit():
            prepend = 'M_'

        #if row['TYPE'] != 'QUAD' :
        #    print name,row['TYPE'],row['K1']

        if row['LENGTH'] == 0 :
            a.AddMarker(prepend + name + '_' + str(eCount))
            continue

        if row['TYPE'] == 'MARK':
            a.AddMarker(prepend + name + '_' + str(eCount))
        ###################################################################################
        elif row['TYPE'] == 'MONI':
            a.AddMarker(prepend + name + '_' + str(eCount))
        ###################################################################################
        elif row['TYPE'] == 'DRIFT':
            a.AddDrift(prepend + name + '_' + str(eCount),
                       length=row['LENGTH'],
                       aper1=defaultBpRadius)
        ###################################################################################
        elif row['TYPE'] == 'QUAD':
            if row['LENGTH'] < 1e-7:
                a.AddMarker(prepend + name + '_' + str(eCount))
            else:
                a.AddQuadrupole(prepend + name + '_' + str(eCount),
                                k1=flip*row['K1']/row['LENGTH'],
                                length=row['LENGTH'],
                                tilt=row['ROTATE'],
                                aper1=defaultBpRadius)
        ###################################################################################
        elif row['TYPE'] == 'BEND' :
            angle = row['ANGLE']
            length = row['LENGTH']
            e1    = 0
            e2    = 0
            fint  = 0
            fintx = 0
            hgap  = 0
            #print name,row['K0'],row['K1'],row['K2']
            kws = {}
            if (e1 != 0):
                kws['e1'] = e1
            if (e1 != 0):
                kws['e2'] = e2
            if (fint != 0):
                kws['fint']  = fint
            if (fintx != 0):
                kws['fintx'] = fintx
            if (hgap != 0):
                kws['hgap'] = hgap

            a.AddDipole(name,'sbend',length=row['LENGTH'],angle=angle,**kws)
        ###################################################################################
        elif row['TYPE'] == 'SEXT':
            l = row['LENGTH']
            if l < 1e-7:
                a.AddDrift(prepend + name + '_' + str(eCount),
                           length=l,
                           aper1=defaultBpRadius)
            else:
                if enableSextupoles:
                    k2in = flip*row['K2']
                else:
                    k2in = 0.0
                a.AddSextupole(prepend + name + '_' + str(eCount),
                               length=l,
                               k2=k2in,
                               aper1=defaultBpRadius)
        else:
            #    print "UNKN> ", row['TYPE']
            try :
                unknDict[row['TYPE']] += 1
            except KeyError :
                unknDict[row['TYPE']] = 1
        nameDict[name] += 1

    print unknDict
    a.AddSampler(samplers)
    a.Write(outputFileName)

    return a


'''
            ###################################################################################
    elif c.type[i] == 'OCTU':
        if c.data[i][c.keys['octu']['l']] > 1e-7:
            a.AddDrift(prepend + c.name[i] + '_' + str(eCount),
                       length=float(c.data[i][c.keys['octu']['l']]),
                       aper1=apertures.aper[i])
        else:
            a.AddMarker(prepend + c.name[i] + '_' + str(eCount))
            ###################################################################################
    elif c.type[i] == 'DECU':
        pass
    ###################################################################################
    elif c.type[i] == 'MULT':
        a.AddMarker(prepend + c.name[i] + '_' + str(eCount))
    ###################################################################################
    elif c.type[i] == 'HKIC':
        a.AddDrift(prepend + c.name[i] + '_' + str(eCount),
                   length=float(c.data[i][c.keys['hkic']['l']]),
                   aper1=float(apertures.aper[i]))
    ###################################################################################
    elif c.type[i] == 'VKIC':
        a.AddDrift(prepend + c.name[i] + '_' + str(eCount),
                   length=float(c.data[i][c.keys['vkic']['l']]),
                   aper1=float(apertures.aper[i]))
    ###################################################################################
    elif c.type[i] == 'SBEN':
        if c.data[i][c.keys['sben']['l']] < 1e-7:
            a.AddMarker(prepend + c.name[i] + '_' + str(eCount))
        else:
            # check for large tilt
            # if enableDipoleTiltTransform and float(c.data[i][c.keys['sben']['tilt']]) > 0.2 :
            #    print "inserting transform 3D"
            #    a.AddTransform3D(c.name[i]+"3dt_in",phi=0,theta=0,psi=float(c.data[i][c.keys['sben']['tilt']]))

            # print "SBEN> ",c.name[i],c.data[i][c.keys['sben']['tilt']]

            # check for poleface
            e1in = 0.0
            e2in = 0.0
            if enableDipolePoleFaceRotation:
                e1in = float(c.data[i][c.keys['sben']['e1']])
                e2in = float(c.data[i][c.keys['sben']['e2']])

            a.AddDipole(prepend + c.name[i] + '_' + str(eCount), 'sbend',
                        length=float(c.data[i][c.keys['sben']['l']]),
                        angle=float(c.data[i][c.keys['sben']['angle']]),
                        aper=float(apertures.aper[i]),
                        e1=e1in,
                        e2=e2in,
                        tilt=float(c.data[i][c.keys['sben']['tilt']]))


            # if enableDipoleTiltTransform and float(c.data[i][c.keys['sben']['tilt']]) > 0.2 :
            #    print "removing transform 3D"
            #    a.AddTransform3D(c.name[i]+"3dt_out",phi=0,theta=0,psi=-float(c.data[i][c.keys['sben']['tilt']]))

            ###################################################################################
    elif c.type[i] == 'LCAV':
        length = float(c.data[i][c.keys['lcav']['l']])
        deltaE = (float(c.data[i][c.keys['lcav']['E']]) - float(c.data[i - 1][c.keys['lcav']['E']])) * 1000  # MeV
        gradient = deltaE / length
        a.AddRFCavity(prepend + c.name[i] + '_' + str(eCount),
                      length=length,
                      gradient=-gradient,
                      aper1=apertures.aper[i])
    ###################################################################################
    elif c.type[i] == 'ECOL':

        if collimator == None:
            # make collimator from mad8 file
            #                print "ECOL> ",c.name[i], "mad8 file"
            print c.data[i][c.keys['ecol']['xsize']], c.data[i][c.keys['ecol']['ysize']]
            if (c.data[i][c.keys['ecol']['xsize']] != 0) and (c.data[i][c.keys['rcol']['ysize']]) != 0:
                a.AddECol(prepend + c.name[i] + '_' + str(eCount),
                          length=float(c.data[i][c.keys['ecol']['l']]),
                          xsize=float(c.data[i][c.keys['ecol']['xsize']]),
                          ysize=float(c.data[i][c.keys['ecol']['ysize']]),
                          material='Copper')
            else:
                a.AddDrift(c.name[i] + '_' + str(eCount), c.data[i][c.keys['ecol']['l']])
        else:
            # make collimator from file
            #                print "ECOL> ",c.name[i], "coll file"
            length = float(c.data[i][c.keys['rcol']['l']])
            xsize = float(collimator.getCollimator(c.name[i])['xsize'])
            ysize = float(collimator.getCollimator(c.name[i])['ysize'])
            mater = collimator.getCollimator(c.name[i])['bdsim_material']
            if xsize != 0 and ysize == 0:
                ysize = 0.2
            if xsize == 0 and ysize != 0:
                xsize = 0.2

            if (xsize != 0 and ysize != 0):
                a.AddECol(prepend + c.name[i] + '_' + str(eCount),
                          length=length,
                          xsize=xsize,
                          ysize=ysize,
                          material=mater,
                          bias=biasList)
            else:
                a.AddDrift(prepend + c.name[i] + '_' + str(eCount), float(c.data[i][c.keys['ecol']['l']]))
                ###################################################################################
    elif c.type[i] == 'RCOL':
        if collimator == None:
            #               print "RCOL> ",c.name[i], "mad8 file"
            print c.data[i][c.keys['rcol']['xsize']], c.data[i][c.keys['rcol']['ysize']]
            if (c.data[i][c.keys['rcol']['xsize']] != 0) and (c.data[i][c.keys['rcol']['ysize']]) != 0:
                a.AddRCol(prepend + c.name[i] + '_' + str(eCount),
                          length=float(c.data[i][c.keys['rcol']['l']]),
                          xsize=float(c.data[i][c.keys['rcol']['xsize']]),
                          ysize=float(c.data[i][c.keys['rcol']['ysize']]),
                          material='Copper')
            else:
                a.AddDrift(prepend + c.name[i] + '_' + str(eCount), c.data[i][c.keys['rcol']['l']])
        else:
            # make collimator from file
            #              print "RCOL> ",c.name[i], "coll file"
            length = float(c.data[i][c.keys['rcol']['l']])
            xsize = float(collimator.getCollimator(c.name[i])['xsize'])
            ysize = float(collimator.getCollimator(c.name[i])['ysize'])
            mater = collimator.getCollimator(c.name[i])['bdsim_material']
            if xsize != 0 and ysize == 0:
                ysize = 0.2
            if xsize == 0 and ysize != 0:
                xsize = 0.2

            if (xsize != 0 and ysize != 0):
                a.AddRCol(prepend + c.name[i] + '_' + str(eCount),
                          length=length,
                          xsize=xsize,
                          ysize=ysize,
                          material=mater,
                          bias=biasList)
            else:
                a.AddDrift(prepend + c.name[i] + '_' + str(eCount), float(c.data[i][c.keys['rcol']['l']]))
                ###################################################################################
'''