# pybdsim.Options - generate BDSIM options
# Version 1.0
# L. Nevay
# laurie.nevay@rhul.ac.uk


def ProtonColliderOptions():
    a = Options()
    a.SetPhysicsList('QGSP_BERT')
    a.SetBeamPipeThickness(5,'mm')
    a.SetOuterDiameter(0.5,'m')
    a.SetTunnelRadius(2,'m')
    a.SetNGenerate(100)

    #Non mandatory options
    a.SetBeamPipeRadius(5,'cm')

    #Tunnel
    a.SetBuildTunnel(True)
    a.SetBuildTunnelFloor(True)
    a.SetTunnelRadius(2,'m')
    a.SetTunnelThickness(50,'cm')
    a.SetTunnelOffsetX(1,'m')
    a.SetTunnelFloorOffset(1,'m')
    return a

def ElectronColliderOptions():
    a = Options()
    a.SetPhysicsList('em')
    a.SetBeamPipeThickness(5,'mm')
    a.SetOuterDiameter(1,'m')
    a.SetTunnelRadius(2,'m')
    a.SetNGenerate(100)

    #Production cuts
    a.SetDefaultRangeCut(0.25,"m")
    a.SetProductionCutElectrons(0.25,"m")
    a.SetProductionCutPositrons(0.25,"m")
    a.SetProductionCutPhotons(0.25,"m")

    #Non mandatory options
    a.SetBeamPipeRadius(5,'cm')

    #Tunnel
    a.SetBuildTunnel(True)
    a.SetBuildTunnelFloor(True)
    a.SetTunnelRadius(2,'m')
    a.SetTunnelThickness(50,'cm')
    a.SetTunnelOffsetX(1,'m')
    a.SetTunnelFloorOffset(1,'m')
    return a

def MinimumStandard():
    a = Options()
    # no specific defaults needed
    return a

class Options(dict):
    def __init__(self,*args,**kwargs):
        dict.__init__(self,*args,**kwargs)

    def SetGeneralOption(self, option, value):
        self[option] = value

    def ReturnOptionsString(self):
        s = ''
        if len(self.keys()) == 0:
            print('No options set - empty string')
            return s
        
        numOptions=0
        for k,v in self.items():
            s += ', \n\t'+str(k)+'='+str(v)
            numOptions += 1
        s += ';'
        s2 = s.split('\n')
        s3 = 'option,\t'+s2[1].replace('\t','').replace('\n','').replace(',','').strip()
        if numOptions == 1:
            s3 += '\n'
        else:
            s3 += ',\n'
        s4 = '\n'.join(s2[2:])
        st = s3+s4
        return st

    def SetNGenerate(self,nparticles=1):
        self['ngenerate'] = nparticles

    def SetPhysicsList(self,physicslist=''):
        physicslistlist = [
            'em',
            'em_low',
            'synchrad',
            'optical',
            'hadronic',
            'hadronichp',
            'qgsp_bert',
            'qgsp_bert_hp',
            'qgsp_bic',
            'qgsp_bic_hp',
            'ftfp_bert',
            'ftfp_bert_hp',
            'decay',
            'muon'
            ]
        if len(physicslist.split()) == 1 :
            if physicslist not in physicslistlist:
                raise ValueError('Unknown physicslist: '+physicslist)
            self['physicsList'] = '"' + str(physicslist) + '"'
        else :
            splitphysicslist = physicslist.split()
            for token in splitphysicslist :
                if token not in physicslistlist :
                    raise ValueError('Unknown physicslist: ' + physicslist)

            self['physicsList'] = '"' + str(physicslist) + '"'

    def SetBeamPipeRadius(self,beampiperadius=5,unitsstring='cm'):
        self['beampipeRadius'] = str(beampiperadius) + '*' +unitsstring

    def SetOuterDiameter(self,outerdiameter=2,unitsstring='m'):
        self['outerDiameter'] = str(outerdiameter) + '*' + unitsstring

    def SetTunnelRadius(self,tunnelradius=2,unitsstring='m'):
        self['tunnelRadius'] = str(tunnelradius) + '*' + unitsstring

    def SetBeamPipeThickness(self,bpt,unitsstring='mm'):
        self['beampipeThickness'] = str(bpt) + '*' + unitsstring

    def SetPipeMaterial(self,bpm):
        self['pipeMaterial'] = '"' + str(bpm) + '"'

    def SetVacuumMaterial(self,vm):
        self['vacMaterial'] = '"' + str(vm) + '"'

    def SetVacuumPressure(self,vp):
        """
        Vacuum pressure in bar
        """
        self['vacuumPressure'] = str(vp)

    def SetBuildTunnel(self,tunnel=False):
        if tunnel == True:
            self['buildTunnel'] = 1
        else:
            self['buildTunnel'] = 0

    def SetBuildTunnelFloor(self,tunnelfloor=False):
        if tunnelfloor == True:
            self['buildTunnelFloor'] = 1
        else:
            self['buildTunnelFloor'] = 0

    def SetTunnelThickness(self,tt=1.0,unitsstring='m'):
        self['tunnelThickness'] = str(tt) + '*' + unitsstring

    def SetSoilThickness(self,st=4.0,unitsstring='m'):
        self['tunnelSoilThickness'] = str(st) + '*' + unitsstring

    def SetTunnelMaterial(self,tm):
        self['tunnelMaterial'] = '"' + str(tm) + '"'

    def SetSoilMaterial(self,sm,):
        self['soilMaterial'] = '"' + str(sm) + '"'

    def SetTunnelOffsetX(self,offset=0.0,unitsstring='m'):
        self['tunnelOffsetX'] = str(offset) + '*' + unitsstring

    def SetTunnelOffsetY(self,offset=0.0,unitsstring='m'):
        self['tunnelOffsetY'] = str(offset) + '*' + unitsstring

    def SetTunnelFloorOffset(self,offset=1.0,unitsstring='m'):
        self['tunnelFloorOffset'] = str(offset) + '*' + unitsstring

    def SetSamplerDiameter(self,radius=10,unitsstring='m'):
        self['samplerDiameter'] = str(radius) + '*' + unitsstring

    def SetBLMRadius(self,radius=5,unitsstring='cm'):
        self['blmRad'] = str(radius) + '*' + unitsstring

    def SetBLMLength(self,length=50,unitsstring='cm'):
        self['blmLength'] = str(length) + '*' + unitsstring

    def SetIncludeIronMagField(self,iron=True):
        if iron == True:
            self['includeIronMagFields'] = 1
        else:
            self['includeIronMagFields'] = 0

    def SetDontSplitSBends(self,dontsplitsbends=False):
        if dontsplitsbends:
            self['dontSplitSBends'] = 1
        else:
            self['dontSplitSBends'] = 0

    def SetDeltaChord(self,dc=0.001,unitsstring='m'):
        self['deltaChord'] = str(dc) + '*' + unitsstring

    def SetDeltaIntersection(self,di=10,unitsstring='nm'):
        self['deltaIntersection'] = str(di) + '*' + unitsstring

    def SetChordStepMinimum(self,csm=1,unitsstring='nm'):
        self['chordStepMinimum'] = str(csm) + '*' + unitsstring

    def SetLengthSafety(self,ls=10,unitsstring='um'):
        self['lengthSafety'] = str(ls) + '*' + unitsstring

    def SetMinimumEpsilonStep(self,mes=10,unitsstring='nm'):
        self['minimumEpsilonStep'] = str(mes) + '*' + unitsstring

    def SetMaximumEpsilonStep(self,mes=1,unitsstring='m'):
        self['maximumEpsilonStep'] = str(mes) + '*' + unitsstring

    def SetDeltaOneStep(self,dos=10,unitsstring='nm'):
        self['deltaOneStep'] = str(dos) + '*' + unitsstring

    def SetMaximumStepLength(self,msl=20,unitsstring='m'):
        self['maximumStepLength'] = str(msl) + '*' + unitsstring

    def SetMaximumTrackingTime(self,mtt=-1,unitsstring='s'):
        self['maximumTrackingTime'] = str(mtt) + '*' + unitsstring

    def SetIntegratorSet(self,integratorSet='"bdsim"'):
        self['integratorSet'] = integratorSet

    def SetThresholdCutCharged(self,tcc=100,unitsstring='MeV'):
        self['thresholdCutCharged'] = str(tcc) + '*' + unitsstring

    def SetThresholdCutPhotons(self,tcp=1,unitsstring='MeV'):
        self['thresholdCutPhotons'] = str(tcp) + '*' + unitsstring

    def SetStopTracks(self,stop=True):
        if stop == True:
            self['stopTracks'] = 1
        else:
            self['stopTracks'] = 0

    def SetStopSecondaries(self,stop=True):
        if stop == True:
            self['stopSecondaries'] = 1
        else:
            self['stopSecondaries'] = 0

    def SetSynchRadiationOn(self,on=True):
        if on == True:
            self['synchRadOn'] = 1
        else:
            self['synchRadOn'] = 0

    def SetTrackSRPhotons(self,track=True):
        if track == True:
            self['srTrackPhotons'] = 1
        else:
            self['srTrackPhotons'] = 0

    def SetSRLowX(self,lowx=True):
        if lowx == True:
            self['srLowX'] = 1
        else:
            self['srLowX'] = 0

    def SetSRMultiplicity(self,srm=2.0):
        self['srMultiplicity'] = srm

    def SetProductionCutPhotons(self,pc=100,unitsstring='keV'):
        self['prodCutPhotons'] = str(pc) + '*' + unitsstring

    def SetProductionCutElectrons(self,pc=100,unitsstring='keV'):
        self['prodCutElectrons'] = str(pc) + '*' + unitsstring

    def SetProductionCutPositrons(self,pc=100,unitsstring='keV'):
        self['prodCutPositrons'] = str(pc) + '*' + unitsstring

    def SetCherenkovOn(self,on=True):
        if on == True:
            self['turnOnCerenkov'] = 1
        else:
            self['tunrOnCerenkov'] = 0

    def SetDefaultRangeCut(self,drc=0.7,unitsstring='mm'):
        self['defaultRangeCut'] = str(drc) + '*' + unitsstring

    def SetGamma2MuonEnahncementFactor(self,ef=2):
        self['gammToMuFe'] = ef

    def SetEPAnnihilation2MuonEnhancementFactor(self,ef=2):
        self['annihiToMuFe'] = ef

    def SetEPAnnihilation2HadronEnhancementFactor(self,ef=2):
        self['eetoHadronsFe'] = ef

    def SetEMLeadParticleBiasing(self,on=True):
        if on == True:
            self['useEMLPB'] = 1
        else:
            self['useEMLPB'] = 0

    def SetLPBFraction(self,fraction=0.5):
        self['LPBFraction'] = fraction

    def SetRandomSeed(self,rs=0):
        self['randomSeed'] = rs

    def SetNGenerate(self,nparticles=10):
        self['ngenerate'] = nparticles

    def SetWritePrimaries(self,on=True):
        if on == True:
            self['writePrimaries'] = 1
        else:
            self['writePrimaries'] = 0

    def SetELossHistBinWidth(self,width):
        self['elossHistoBinWidth'] = width

    def SetSensitiveBeamlineComponents(self,on=True):
        if True:
            self['sensitiveBeamLineComponents'] = 1
        else:
            self['sensitiveBeamLineComponents'] = 0

    def SetSensitiveBeamPipe(self,on=True):
        if True:
            self['sensitiveBeamPipe'] = 1
        else:
            self['sensitiveBeamPipe'] = 0

    def SetSenssitiveBLMs(self,on=True):
        if True:
            self['sensitiveBLMs'] = 1
        else:
            self['sensitiveBLMs'] = 0

    def SetStoreTrajectory(self,on=True):
        if True:
            self['storeTrajectory'] = 1
        else:
            self['storeTrajectory'] = 0

    def SetStoreTrajectoryParticle(self, particle = "muon"):
        self['storeTrajectoryParticle'] = particle

    def SetMagnetGeometryType(self, magnetGeometryType = '"none"'):
        self['magnetGeometryType'] = magnetGeometryType

    def SetStoreMuonTrajectory(self,on=True):
        if True:
            self['storeMuonTrajectory'] = 1
        else:
            self['storeMuonTrajectory'] = 0

    def SetStoreNeutronTrajectory(self,on=True):
        if True:
            self['storeNeutronTrajectory'] = 1
        else:
            self['storeNeutronTrajectory'] = 0

    def SetTrajectoryCutGTZ(self,gtz=0.0,unitsstring='m'):
        self['trajCutGTZ'] = str(gtz) + '*' + unitsstring

    def SetTrajectoryCutLTR(self,ltr=10.0,unitsstring='m'):
        self['trajCutLTR'] = str(ltr) + '*' + unitsstring

    def SetPrintModuloFraction(self,pmf=1e-2):
        self['printModuloFraction'] = pmf

    def SetNPerFile(self,nperfile=100):
        self['nperfile'] = nperfile

    def SetNLinesIgnore(self,nlines=0):
        self['nlinesIgnore'] = nlines

    def SetIncludeFringeFields(self,on=True):
        if on == True:
            self['includeFringeFields'] = 1
        else:
            self['includeFringeFields'] = 0

    def SetDefaultBiasVaccum(self, biases=""):
        self["defaultBiasVacuum"] = biases

    def SetDefaultBiasMaterial(self, biases=""):
        self["defaultBiasMaterial"] = biases

class Editor :
    def __init__(self, fileName) :
        self.fileName = fileName
