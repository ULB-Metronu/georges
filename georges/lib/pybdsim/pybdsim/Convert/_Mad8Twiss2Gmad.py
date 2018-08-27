#!/usr/bin/env python2.6

import numpy as _np
import optparse as _op
from collections import OrderedDict

import pymad8

from .. import Builder 
from .. import Beam
from pybdsim.Options import Options
import pybdsim.XSecBias


def Mad8Twiss2Gmad(inputFileName, outputFileName, 
                   istart                       = 0,
                   beam                         = ["nominal"],
                   gemit                        = (1e-10,1e-10), 
                   mad8FileName                 = "",                  
                   collimator                   = "collimator.dat", 
                   apertures                    = "apertures.dat",
                   samplers                     = 'all',
                   options                      = True,
                   flip                         = 1, 
                   enableSextupoles             = True, 
                   enableOctupoles              = True,
                   enableDecapoles              = True,
                   openApertures                = True,
                   openCollimators              = True,
                   enableDipoleTiltTransform    = True,
                   enableDipolePoleFaceRotation = True,
                   enableSr                     = False,
                   enableSrScaling              = False,
                   enableMuon                   = False,
                   enableMuonBias               = True):
    """
    Convert MAD8 twiss output to a BDSIM model in GMAD syntax.

    """

    # open mad output
    o = pymad8.Output.OutputReader()
    c, t = o.readFile(inputFileName,'twiss')

    # check type of start 
    # if string find element number
    if type(istart) == str : 
        print "Mad8Twiss2Gmad> finding start : ",istart 
        istart = c.findByName(istart)[0]
        print "Mad8Twiss2Gmad> using index   : ",istart
        
    # load Collimator db or use instance
    if type(collimator) == str : 
        collimator = Mad8CollimatorDatabase(collimator) 
        if openCollimators :
            collimator.openCollimators()
    # load Aperture db or use instance
    if type(apertures) == str : 
        apertures = Mad8ApertureDatabase(apertures) 
        if openApertures :
            apertures.openApertures()

    print collimator

    # Need nominal energy for SR calculations
    energy = c.data[istart][c.keys['drif']['E']]

    # create machine instance
    # TODO : Need to extract nominal energy from file
    a = Builder.Machine(sr=enableSrScaling, energy0=float(energy))
    
    # load mad8
    if mad8FileName != "" : 
        particle = 'e+'
        m8 = pymad8.Output.Mad8(mad8FileName)
        if m8.particle == 'ELECTRON' :
            particle = 'e-'
            flip     = -1

        elif m8.particle == 'POSITRON' : 
            particle = 'e+'
            flip     = 1
    else : 
        particle = 'e-'
        flip     = -1

    # create beam (emit and energy spread)
    esprd = 0.0
    if type(gemit) == str : 
        echoVals = pymad8.Output.EchoValue(gemit)
        echoVals.loadValues()
        gemit = _np.zeros(2)
        gemit[0] = echoVals.valueDict['EMITX']
        gemit[1] = echoVals.valueDict['EMITY']
        esprd    = echoVals.valueDict['ESPRD']

    # create beam
    beamname = beam[0]
    if beamname == "reference" :
        b = Beam.Beam(particle, energy, "reference")
        a.AddBeam(b)
    elif beamname == "nominal" :
        b = Mad8Twiss2Beam(t,istart,particle,energy)
        b._SetEmittanceX(gemit[0],'m')
        b._SetEmittanceY(gemit[1],'m')
        b._SetSigmaE(esprd)
        a.AddBeam(b)
    elif beamname == "halo" :
        b = Mad8Twiss2Beam(t,istart,particle,energy)
        b.SetDistributionType("halo")
        betx = t.data[istart][t.keys['betx']]
        bety = t.data[istart][t.keys['bety']]
        alfx = t.data[istart][t.keys['alfx']]
        alfy = t.data[istart][t.keys['alfy']]

        # 5 13
        # 36 92
        nsigxmin = beam[1]
        nsigxmax = beam[2]
        nsigymin = beam[3]
        nsigymax = beam[4]

        b._SetHaloEmittanceX(nsigxmin**2*gemit[0],'m')
        b._SetHaloEmittanceY(nsigymin**2*gemit[1],'m')
        b._SetHaloEnvelopeEmitX(nsigxmax**2*gemit[0],'m')
        b._SetHaloEnvelopeEmitY(nsigymax**2*gemit[1],'m')
        b._SetEnvelopeX(_np.sqrt(nsigxmax**2*betx*gemit[0]))
        b._SetEnvelopeXP(_np.sqrt(nsigxmax**2*(1+alfx**2)/betx*gemit[0]))
        b._SetEnvelopeY(_np.sqrt(nsigymax**2*betx*gemit[1]))
        b._SetEnvelopeYP(_np.sqrt(nsigymax**2*(1+alfy**2)/bety*gemit[1]))
        b._SetHaloPSWeightParameter(-1.0)
        b._SetHaloPSWeightFunction("oneoverr")
        b._SetSigmaE(esprd)
        a.AddBeam(b)

    # create options 
    if options : 
        o = Options.ElectronColliderOptions()
        o.SetBuildTunnel(False)
        o.SetBuildTunnelFloor(False)
        o.SetMagnetGeometryType("none")
        o.SetPrintModuloFraction(1e-5)

        process = 'em'
        if enableSr :
            process += ' synchrad'
        if enableMuon :
            process += ' muon'
            o.SetStoreTrajectory(True)
            o.SetStoreTrajectoryParticle("mu+ mu-")


        o.SetPhysicsList(process)
        a.AddOptions(o)

    # create bias options
    biasList = ""
    if enableMuonBias :
        gmb = XSecBias.XSecBias("gmb","gamma","GammaToMuPair","1e2","1")
        pmb = XSecBias.XSecBias("pmb","e+","AnnihiToMuPair","1e2","1")
        a.AddBias(gmb)
        a.AddBias(pmb)

        biasList = "gmb pmb"



    # iterate through objects and build machine
    for i in range(istart,len(c.name),1) : 
        # unique(c.type)
        # ['', 'BLMO', 'DRIF', 'ECOL', 'HKIC', 'IMON', 'INST', 'LCAV', 'MARK',
        # 'MATR', 'MONI', 'PROF', 'QUAD', 'RCOL', 'SBEN', 'SOLE', 'VKIC',
        #       'WIRE']

        # print element
        # print i,c.name[i],c.type[i]

        # check name dictionary 
        try : 
            eCount = nameDict[c.name[i]]
        except KeyError : 
            nameDict[c.name[i]] = 0
            eCount = nameDict[c.name[i]]

        # check if name starts with a number 
        prepend = ''
        if c.name[i][0].isdigit() :
            prepend = 'M_'

#        print c.name[i]+'_'+str(eCount)

        if c.type[i] == '' : 
            a.AddMarker(prepend+c.name[i]+'_'+str(eCount))
###################################################################################
        elif c.type[i] == 'DRIF' : 
            a.AddDrift(prepend+c.name[i]+'_'+str(eCount),
                       length=float(c.data[i][c.keys['drif']['l']]), 
                       aper1=float(apertures.aper[i]))
###################################################################################
        elif c.type[i] == 'MARK' : 
            a.AddMarker(prepend+c.name[i]+'_'+str(eCount)) 
###################################################################################
        elif c.type[i] == 'SOLE' : 
            a.AddMarker(prepend+c.name[i]+'_'+str(eCount))
###################################################################################
        elif c.type[i] == 'INST' : 
            a.AddMarker(prepend+c.name[i]+'_'+str(eCount)) 
###################################################################################
        elif c.type[i] == 'MONI' : 
            a.AddMarker(prepend+c.name[i]+'_'+str(eCount))
###################################################################################
        elif c.type[i] == 'IMON' :
            a.AddMarker(prepend+c.name[i]+'_'+str(eCount))
###################################################################################
        elif c.type[i] == 'BLMO' : 
            a.AddMarker(prepend+c.name[i]+'_'+str(eCount))
###################################################################################
        elif c.type[i] == 'WIRE' :
            a.AddMarker(prepend+c.name[i]+'_'+str(eCount))
###################################################################################
        elif c.type[i] == 'QUAD' : 
            if c.data[i][c.keys['quad']['l']] < 1e-7 : 
                a.AddMarker(prepend+c.name[i]+'_'+str(eCount))
            else : 
                a.AddQuadrupole(prepend+c.name[i]+'_'+str(eCount),
                                k1     = float(c.data[i][c.keys['quad']['k1']])*flip,
                                length = float(c.data[i][c.keys['quad']['l']]),
                                tilt   = float(c.data[i][c.keys['quad']['tilt']]),
                                aper1  = float(apertures.aper[i]))
###################################################################################
        elif c.type[i] == 'SEXT' : 
            l = float(c.data[i][c.keys['sext']['l']])
            if l < 1e-7 :
                a.AddDrift(prepend+c.name[i]+'_'+str(eCount),
                           length=l,
                           aper1=apertures.aper[i])
            else : 
                if enableSextupoles : 
                    k2in=float(c.data[i][c.keys['sext']['k2']])*flip
                else : 
                    k2in=0.0                    
                a.AddSextupole(prepend+c.name[i]+'_'+str(eCount),
                               length=l,
                               k2=k2in,
                               aper1=apertures.aper[i])
###################################################################################
        elif c.type[i] == 'OCTU' : 
            if c.data[i][c.keys['octu']['l']] > 1e-7 : 
                a.AddDrift(prepend+c.name[i]+'_'+str(eCount),
                           length=float(c.data[i][c.keys['octu']['l']]),
                           aper1=apertures.aper[i])
            else : 
                a.AddMarker(prepend+c.name[i]+'_'+str(eCount))
###################################################################################
        elif c.type[i] == 'DECU' : 
            pass
###################################################################################
        elif c.type[i] == 'MULT' : 
                a.AddMarker(prepend+c.name[i]+'_'+str(eCount))
###################################################################################
        elif c.type[i] == 'HKIC' : 
            a.AddDrift(prepend+c.name[i]+'_'+str(eCount),
                       length=float(c.data[i][c.keys['hkic']['l']]),
                       aper1=float(apertures.aper[i]))
###################################################################################
        elif c.type[i] == 'VKIC' : 
            a.AddDrift(prepend+c.name[i]+'_'+str(eCount),
                       length=float(c.data[i][c.keys['vkic']['l']]),
                       aper1=float(apertures.aper[i]))
###################################################################################
        elif c.type[i] == 'SBEN' : 
            if c.data[i][c.keys['sben']['l']] < 1e-7 : 
                a.AddMarker(prepend+c.name[i]+'_'+str(eCount))
            else : 
                # check for large tilt                
                #if enableDipoleTiltTransform and float(c.data[i][c.keys['sben']['tilt']]) > 0.2 :
                #    print "inserting transform 3D"
                #    a.AddTransform3D(c.name[i]+"3dt_in",phi=0,theta=0,psi=float(c.data[i][c.keys['sben']['tilt']]))
                    
                # print "SBEN> ",c.name[i],c.data[i][c.keys['sben']['tilt']]

                # check for poleface
                e1in = 0.0
                e2in = 0.0
                if enableDipolePoleFaceRotation :
                    e1in = float(c.data[i][c.keys['sben']['e1']])
                    e2in = float(c.data[i][c.keys['sben']['e2']])

                a.AddDipole(prepend+c.name[i]+'_'+str(eCount),'sbend',
                            length= float(c.data[i][c.keys['sben']['l']]), 
                            angle = float(c.data[i][c.keys['sben']['angle']]),
                            aper  = float(apertures.aper[i]),
                            e1    = e1in,
                            e2    = e2in,
                            tilt  = float(c.data[i][c.keys['sben']['tilt']]))


                #if enableDipoleTiltTransform and float(c.data[i][c.keys['sben']['tilt']]) > 0.2 :
                #    print "removing transform 3D"
                #    a.AddTransform3D(c.name[i]+"3dt_out",phi=0,theta=0,psi=-float(c.data[i][c.keys['sben']['tilt']]))

###################################################################################
        elif c.type[i] == 'LCAV' : 
            length   = float(c.data[i][c.keys['lcav']['l']])
            deltaE   = (float(c.data[i][c.keys['lcav']['E']])-float(c.data[i-1][c.keys['lcav']['E']]))*1000 # MeV 
            gradient = deltaE/length
            a.AddRFCavity(prepend+c.name[i]+'_'+str(eCount),
                          length=length, 
                          gradient=-gradient,
                          aper1=apertures.aper[i])
###################################################################################
        elif c.type[i] == 'ECOL' :

            if collimator == None : 
                # make collimator from mad8 file
#                print "ECOL> ",c.name[i], "mad8 file"
                print c.data[i][c.keys['ecol']['xsize']],c.data[i][c.keys['ecol']['ysize']]
                if (c.data[i][c.keys['ecol']['xsize']] != 0) and (c.data[i][c.keys['rcol']['ysize']]) != 0 : 
                    a.AddECol(prepend+c.name[i]+'_'+str(eCount), 
                              length  = float(c.data[i][c.keys['ecol']['l']]), 
                              xsize   = float(c.data[i][c.keys['ecol']['xsize']]), 
                              ysize   = float(c.data[i][c.keys['ecol']['ysize']]),
                              material= 'Copper')
                else : 
                    a.AddDrift(c.name[i]+'_'+str(eCount),c.data[i][c.keys['ecol']['l']])
            else : 
                # make collimator from file
#                print "ECOL> ",c.name[i], "coll file"
                length=float(c.data[i][c.keys['rcol']['l']])
                xsize = float(collimator.getCollimator(c.name[i])['xsize'])
                ysize = float(collimator.getCollimator(c.name[i])['ysize'])
                mater = collimator.getCollimator(c.name[i])['bdsim_material']
                if xsize != 0 and ysize == 0 :
                    ysize = 0.2
                if xsize == 0 and ysize != 0 :
                    xsize = 0.2

                if (xsize != 0 and ysize != 0):
                    a.AddECol(prepend+c.name[i]+'_'+str(eCount), 
                              length  = length,
                              xsize   = xsize,
                              ysize   = ysize,
                              material= mater,
                              bias    = biasList)
                else : 
                    a.AddDrift(prepend+c.name[i]+'_'+str(eCount),float(c.data[i][c.keys['ecol']['l']]))
###################################################################################
        elif c.type[i] == 'RCOL' :
            if collimator == None : 
#               print "RCOL> ",c.name[i], "mad8 file"
                print c.data[i][c.keys['rcol']['xsize']],c.data[i][c.keys['rcol']['ysize']]
                if (c.data[i][c.keys['rcol']['xsize']] != 0) and (c.data[i][c.keys['rcol']['ysize']]) != 0 : 
                    a.AddRCol(prepend+c.name[i]+'_'+str(eCount), 
                              length  = float(c.data[i][c.keys['rcol']['l']]), 
                              xsize   = float(c.data[i][c.keys['rcol']['xsize']]), 
                              ysize   = float(c.data[i][c.keys['rcol']['ysize']]),
                              material= 'Copper')            
                else : 
                    a.AddDrift(prepend+c.name[i]+'_'+str(eCount),c.data[i][c.keys['rcol']['l']])
            else : 
                # make collimator from file
 #              print "RCOL> ",c.name[i], "coll file"
                length= float(c.data[i][c.keys['rcol']['l']])
                xsize = float(collimator.getCollimator(c.name[i])['xsize'])
                ysize = float(collimator.getCollimator(c.name[i])['ysize'])
                mater = collimator.getCollimator(c.name[i])['bdsim_material']
                if xsize != 0 and ysize == 0 :
                    ysize = 0.2
                if xsize == 0 and ysize != 0 :
                    xsize = 0.2

                if (xsize != 0 and ysize != 0) :
                    a.AddRCol(prepend+c.name[i]+'_'+str(eCount),
                              length  = length,
                              xsize   = xsize,
                              ysize   = ysize,
                              material= mater,
                              bias=biasList)
                else : 
                    a.AddDrift(prepend+c.name[i]+'_'+str(eCount),float(c.data[i][c.keys['rcol']['l']]))
###################################################################################
        else :
            print "UNKN> ",c.type[i]
###################################################################################

        nameDict[c.name[i]] += 1 

    a.AddSampler(samplers)
    a.Write(outputFileName)

    return a

def Mad8Twiss2Beam(t, istart, particle, energy) :            
    
    betx = t.data[istart][t.keys['betx']]
    bety = t.data[istart][t.keys['bety']]
    alfx = t.data[istart][t.keys['alfx']]
    alfy = t.data[istart][t.keys['alfy']]

    beam = Beam.Beam(particle,energy,'gausstwiss')
    beam._SetBetaX(betx)
    beam._SetBetaY(bety)
    beam._SetAlphaX(alfx)
    beam._SetAlphaY(alfy)

    return beam


def Mad8MakeOptions(inputTwissFile, inputEchoFile) :
    # open mad output
    o = pymad8.Output.OutputReader()
    c, t = o.readFile(inputFileName,'twiss')
    
    # get initial beam pipe size
    a = c.getApertures(raw=False)
#    a = c.getApertures(raw=True)

    # get values from echo of mad8 output (particle type, beam energy, energy spread)
    echoVals = pymad8.Output.EchoValue(inputEchoFile)
    echoVals.loadValues()

def Mad8MakeApertureTemplate(inputFileName, outputFileName="apertures_template.dat") : 
    # open mad output
    o = pymad8.Output.OutputReader()
    c, t = o.readFile(inputFileName,'twiss')
    a = c.getApertures(raw=False)
#    a = c.getApertures(raw=True)

    # write apertures to file
    f = open(outputFileName,"w")     

    for i in range(0,len(c.name),1) : 
        f.write(c.name[i]+' '+str(a[i])+'\n');
    f.close()

def Mad8MakeCollimatorTemplate(inputFileName,outputFileName="collimator_template.dat") : 
    '''
    Read Twiss file and generate template of collimator file 
    inputFileName  = "twiss.tape"
    outputFileName = "collimator.dat"
    collimator.dat must be edited to provide types and materials.py, apertures will be defined from lattice
    '''
    # open mad output
    o = pymad8.Output.OutputReader()
    c, t = o.readFile(inputFileName,'twiss')

    # open collimator output file
    f = open(outputFileName,"w") 

    for i in range(0,len(c.name),1) : 
        if c.type[i] == 'RCOL' :
            f.write(c.name[i]+"\t"+"TYPE"+"\t"+str(c.data[i][c.keys['rcol']['l']])+"\t"+str(c.data[i][c.keys['rcol']['xsize']])+"\t"+str(c.data[i][c.keys['rcol']['ysize']])+"\t"+"Copper"+"\t"+"GEOM"+"\t"+"SIGMA"+"\n")
        if c.type[i] == 'ECOL' : 
            f.write(c.name[i]+"\t"+"TYPE"+"\t"+str(c.data[i][c.keys['ecol']['l']])+"\t"+str(c.data[i][c.keys['ecol']['xsize']])+"\t"+str(c.data[i][c.keys['ecol']['ysize']])+"\t"+"Copper"+"\t"+"GEOM"+"\t"+"SIGMA"+"\n")
    f.close()
        
class Mad8ApertureDatabase: 
    def __init__(self,apertureFileName) : 
        self.apertureFileName = apertureFileName
        self.loadApertures(self.apertureFileName) 

    def loadApertures(self,fileName) :
        f = open(fileName)
        self.name = [] 
        self.aper = [] 
        for l in f : 
            t = l.split() 
            self.name.append(t[0])
            self.aper.append(float(t[1]))

    def openApertures(self,size = 0.1) : 
        for i in range(0,len(self.aper)):
            self.aper[i] = size
        
class Mad8CollimatorDatabase: 
    '''
    Load collimator file into memory and functions to open and manipulate collimator system
    c = Mad8CollimatorDataBase(fileName)
    fileName = "collimator.dat"
    file format
    <name> <type> <length> <x_size/2> <ysize/2> <material> <geom>
    <length> includes the tapers, so entire length along machine 
    '''

    def __init__(self,collimatorFileName) :
        self.collimatorFileName = collimatorFileName
        self.loadCollimatorDb(self.collimatorFileName) 
        
    def loadCollimatorDb(self,collimatorFileName) : 
        f = open(collimatorFileName) 

        inx = 0 

        self._coll = OrderedDict()
        self._collNames = []
        for l in f : 
            t = l.split()
            name     = t[0]            
            type     = t[1]
            length   = float(t[2])
            xsize    = float(t[3])
            ysize    = float(t[4]) 
            material = t[5]
            geom     = t[6]
            setting = t[7]
            inx = inx + 1 

            d = {'index':inx, 'type':type, 'l':length, 'xsize':xsize,
                 'ysize':ysize, 'bdsim_material':material, 'bdsim_geom':geom, 'setting':setting}

            self._coll[name] = d     
            self._collNames.append(name)

    def openCollimators(self,openHalfSizeX=0.2, openHalfSizeY=0.2) :
        for c in self._coll.keys() :
            self._coll[c]['xsize'] = openHalfSizeX
            self._coll[c]['ysize'] = openHalfSizeY

    def setCollimator(self,collimator,halfSizeX,halfSizeY) : 
        self._coll[collimator]['xsize'] = halfSizeX
        self._coll[collimator]['ysize'] = halfSizeY
    
    def getCollimators(self) : 
#        return self._coll.keys()
        return self._collNames

    def getCollimator(self, name) : 
        return self._coll[name]

    def getDict(self) : 
        return self._coll

    def __str__(self) : 
        s = ''
        for k in self.getCollimators() : 
            s += k+"\t"+self._coll[k]['type']+"\t"+str(self._coll[k]['l'])+"\t"+str(self._coll[k]['xsize'])+"\t"+str(self._coll[k]['ysize'])+"\t"+self._coll[k]['bdsim_material']+"\t"+self._coll[k]['bdsim_geom']+"\t"+self._coll[k]['setting']+"\n"
        return s

    def write(self,fileName) : 
        f = open(fileName,"w")
        
        for k in self.getCollimators() : 
            f.write(k+"\t"+self._coll[k]['type']+"\t"+str(self._coll[k]['l'])+"\t"+str(self._coll[k]['xsize'])+"\t"+str(self._coll[k]['ysize'])+"\t"+self._coll[k]['bdsim_material']+"\t"+self._coll[k]['bdsim_geom']+"\t"+self._coll[k]['setting']+"\n")
        

def main() : 
    usage = "usage : %prog [inputFileName]" 
    parser = _op.OptionParser(usage)
    parser.add_option("-i","--input",  action="store",   type="string",     dest="inputFileName",                           help="Input file name")
    parser.add_option("-o","--ouput",  action="store",   type="string",     dest="outputFileName",default="output",         help="output base file name")
    parser.add_option("-s","--start",  action="store",   type="string",     dest="start",         default="0",              help="starting element (named or index)")
    parser.add_option("-b","--beam",   action="store_true",                 dest="beam",          default=True,             help="generate beam") 
    parser.add_option("-x","--emitx",  action="store",   type="float",      dest="gemitx",        default=1e-10,            help="geometric emittance in x")
    parser.add_option("-y","--emity",  action="store",   type="float",      dest="gemity",        default=1e-10,            help="geometric emittance in x")
    parser.add_option("-c","--coll",   action="store",   type="string",     dest="coll",          default="collimator.dat", help="collimator defn file")
    parser.add_option("-a","--sampler",action="store",   type="string",     dest="samplers",      default="all",            help="samplers (all|)") 

    (options, args) = parser.parse_args()

    if options.inputFileName == None : 
        print "_Mad8Twiss2Gmad> Require input file name"
        parser.print_usage()
        return 
    print '_Mad8Twiss2Gmad> inputFileName,outputFileName,start,samplers,beam,gemitx,gemity'
    print '_Mad8Twiss2Gmad>', options.inputFileName,options.outputFileName,options.start,options.samplers,options.beam,options.gemitx,options.gemity
    
    # try to decode the start point either value or name
    try :
        options.start = int(options.start) 
    except ValueError : 
        pass 
    
    Mad8Twiss2Gmad(options.inputFileName, options.outputFileName, options.start, options.beam, (options.gemitx,options.gemity),options.coll,options.samplers)
    
if __name__ == "__main__":
    main()
