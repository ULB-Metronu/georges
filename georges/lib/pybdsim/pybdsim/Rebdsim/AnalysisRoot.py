import ROOT             as _ROOT
from ROOT import TChain as _TChain
from ROOT import TH1F   as _TH1F
from ROOT import TFile  as _TFile

import matplotlib.pyplot as _plt

class AnalysisRoot:
    def __init__(self,filelist):
        if filelist != None:
            self.Load(filelist)

    def Load(self,filelist):
        self.filelist = filelist
        self.nfiles   = 0
        self._LoadPrimaries()
        self._LoadElossTree()
        self._LoadPlossTree()
        self._LoadSamplerTree()
        if self.nfiles == 0:
            self._NoFilesWarning()

    def _LoadPrimaries(self) : 
        self.primaries = _TChain("primaries") 
        if type(self.filelist) == list : 
            for f in self.filelist : 
                self.primaries.Add(f)
                
    def _LoadSamplerTree(self) : 
        # open a single file and determine the available samplers 
        f = _TFile.Open(self.filelist[0])
        k = f.GetListOfKeys()

        self.samplerNames = [] 
        self.samplerDict  = {}        

        for i in  range(0,k.GetEntries(),1) : 
            if k[i].GetName().find("Sampler") == 0 : 
                self.samplerNames.append(k[i].GetName()) 

        # Sampler chain dictionary
        for n in self.samplerNames : 
            self.samplerDict[n] = _TChain(n)
            for f in self.filelist : 
                self.samplerDict[n].Add(f)

    def _LoadElossTree(self):
        self.elossch = _TChain("ElossTree")
        if type(self.filelist) == list:
            for f in self.filelist:
                self.elossch.Add(f)
                self.nfiles += 1
        else:
            self.nfiles = self.elossch.Add(self.filelist)
        
    def _LoadPlossTree(self):
        self.plossch = _TChain("PlossTree")
        if type(self.filelist) == list:
            for f in self.filelist:
                self.plossch.Add(f)
                self.nfiles += 1
        else:
            self.nfiles = self.plossch.Add(self.filelist)

    def _NoFilesWarning(self):
        raise UserWarning("No root files were loaded")
        
    def PlotEnergyLossHistogram(self,nbins=100,xmin=0,xmax=100,normalised=False,weighted=True):
        self.elossch.SetBranchStatus("s",1)
        if weighted:
            self.elossch.SetBranchStatus("E",1)
        self.elosshist = _TH1F("hist","Energy Loss",nbins,xmin,xmax)

        #fill histogram
        nentries = int(self.elossch.GetEntries())
        for i in xrange(nentries):
            self.elossch.GetEntry(i)
            if weighted:
                self.elosshist.Fill(self.elossch.s,self.elossch.E)
            else:
                self.elosshist.Fill(self.elossch.s)
        
        #normalise if required
        if normalised:
            scale = 1./max(self.elosshist)
            print 'Normalising by factor',scale
            self.normalisation = scale
            self.elosshist.Scale(scale)
        
        _PlotHistogram(self.elosshist,normalised)

    def PlotPrimaryLossHistogram(self,nbins=100,xmin=0,xmax=100,normalised=True,weighted=False):
        self.plossch.SetBranchStatus("s",1)
        if weighted:
            self.plossch.SetBranchStatus("E",1)
        self.plosshist = _TH1F("hist","Primary Losses",nbins,xmin,xmax)

        #fill histogram
        nentries = int(self.plossch.GetEntries())
        for i in xrange(nentries):
            self.plossch.GetEntry(i)
            if weighted:
                self.plosshist.Fill(self.plossch.s,self.plossch.E)
            else:
                self.plosshist.Fill(self.plossch.s)
        
        #normalise if required
        if normalised:
            scale = 1./max(self.plosshist)
            self.normalisation = scale
            print 'Normalising by factor',scale
            self.plosshist.Scale(scale)
        
        _PlotHistogram(self.plosshist,normalised)
        

def _PlotHistogram(hist,normalised,title=''):
    #plot histogram
    fig = _plt.figure(figsize=(12,5))
    ax  = fig.add_subplot(111)
    ax.plot(hist)
    #set x limits to add 2% on either side
    xr = len(hist)
    ax.set_xlim(0-0.02*xr,xr+0.02*xr)
    #set large log yscale
    if not normalised:
        ax.set_ylim(min(hist)*0.7,max(hist*1.4))
    else:
        ax.set_ylim(0.7e-7,1.4)
    
    ax.set_xlabel('S Position (m)')
    if not normalised:
        ax.set_ylabel('Number')
    else:
        ax.set_ylabel('Fraction')
    ax.set_yscale('log')
    ax.grid()
    ax.set_title(title)
    _plt.subplots_adjust(left=0.08,right=0.97)