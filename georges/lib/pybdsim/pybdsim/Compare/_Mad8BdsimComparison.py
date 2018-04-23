import pickle as _pkl
import pylab  as _pl
import pymad8 as _pymad8
import pybdsim as _pybdsim
import matplotlib.pyplot as _plt
import numpy as _np

class Mad8Bdsim :
    def __init__(self, 
                 bdsimFileName = "output.pickle",
                 mad8TwissFileName   = "ebds1.twiss",
                 mad8EnvelopeFileName = "ebds1.envelope") : 
        # load bdsim data
        if bdsimFileName.find("pickle") != -1 :
            f = open(bdsimFileName) 
            self.bdsimData = _pkl.load(f)
            self.bdsimOptics = self.bdsimData['optics']
            f.close()
        elif bdsimFileName.find(".root") != -1 :
            import ROOT as _ROOT 
            import root_numpy as _root_numpy
            f = _ROOT.TFile(bdsimFileName)
            t = f.Get("optics")
            self.bdsimOptics = _root_numpy.tree2rec(t)
            


        # load mad8 data
        r = _pymad8.Mad8.OutputReader()    
        [self.mad8Comm,self.mad8Envel] = r.readFile(mad8EnvelopeFileName,"envel")
        [self.mad8Comm,self.mad8Twiss] = r.readFile(mad8TwissFileName,"twiss")

    def plotSigma(self) : 
        figure = _plt.figure(figsize=(11.6, 7.2))
        gs  = _plt.GridSpec(3,1,height_ratios=[1,3,3])
        ax0 = figure.add_subplot(gs[0],projection="_My_Axes")
        _pymad8.Plot.drawMachineLattice(self.mad8Comm,self.mad8Twiss)

        ax1 = _plt.subplot(gs[1])
        _pl.plot(self.mad8Envel.getColumn('suml'),_pl.sqrt(self.mad8Envel.getColumn('s11'))*1e6,"+-",label="MAD8")
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Sigma_x']/1e-6,"x--",label="BDSIM")
        #_pl.xlim(0,2275)
        #_pl.ylim(0,3000)
        _pl.legend(loc=0)
        _pl.ylabel("$\\sigma_x$ [$\mu$m]")

        ax2 = _plt.subplot(gs[2])
        _pl.plot(self.mad8Envel.getColumn('suml'),_pl.sqrt(self.mad8Envel.getColumn('s33'))*1e6,"+-")
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Sigma_y']/1e-6,"x--")
        #_pl.xlim(0,2275)
        #_pl.ylim(0,100)
        _pl.ylabel("$\\sigma_y$ [$\mu$m]")
        _pl.xlabel("$S$ [m]")

        _pymad8.Plot.setCallbacks(figure,ax0,[ax1,ax2],self.mad8Twiss)

        _pl.savefig("mad8bdsim_sigma.pdf")

        
    def plotSigmaPrim(self) : 
        figure = _plt.figure(figsize=(11.6, 7.2))
        gs  = _plt.GridSpec(3,1,height_ratios=[1,3,3])
        ax0 = figure.add_subplot(gs[0],projection="_My_Axes")
        _pymad8.Plot.drawMachineLattice(self.mad8Comm,self.mad8Twiss)

        ax1 = _plt.subplot(gs[1])
        _pl.plot(self.mad8Envel.getColumn('suml'),_pl.sqrt(self.mad8Envel.getColumn('s22'))*1e6,"+-",label="MAD8")
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Sigma_xp']/1e-6,"x--",label="BDSIM")
        #_pl.xlim(0,2275)
        #_pl.ylim(0,3000)
        _pl.legend(loc=0)
        _pl.ylabel("$\\sigma^{'}_{x}$ [$\mu$m]")

        ax2 = _plt.subplot(gs[2])
        _pl.plot(self.mad8Envel.getColumn('suml'),_pl.sqrt(self.mad8Envel.getColumn('s44'))*1e6,"+-")
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Sigma_yp']/1e-6,"x--")
        #_pl.xlim(0,2275)
        #_pl.ylim(0,100)
        _pl.ylabel("$\\sigma^{'}_{y}$ [$\mu$m]")
        _pl.xlabel("$S$ [m]")

        _pymad8.Plot.setCallbacks(figure,ax0,[ax1,ax2],self.mad8Twiss)

        _pl.savefig("mad8bdsim_sigma_prim.pdf")

    def plotOrbit(self) :
        figure = _plt.figure(figsize=(11.6, 7.2))
        gs  = _plt.GridSpec(3,1,height_ratios=[1,3,3])
        ax0 = figure.add_subplot(gs[0],projection="_My_Axes")
        _pymad8.Plot.drawMachineLattice(self.mad8Comm,self.mad8Twiss)

        ax1 = _plt.subplot(gs[1])
        _pl.plot(self.mad8Envel.getColumn('suml'), _np.zeros(len(self.mad8Envel.getColumn('suml'))),"+-",label="MAD8") #mad8 orbit perfectly on reference
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Mean_x']/1e-6,"x--",label="BDSIM")
        #_pl.xlim(0,2275)
        #_pl.ylim(0,3000)
        _pl.legend(loc=0)
        _pl.ylabel("$\\overline{x}$ [$\mu$m]")

        ax2 = _plt.subplot(gs[2])
        _pl.plot(self.mad8Envel.getColumn('suml'), _np.zeros(len(self.mad8Envel.getColumn('suml'))) ,"+-")
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Mean_y']/1e-6,"x--")
        _pl.xlim(0,2275)
        #_pl.ylim(0,100)
        _pl.ylabel("$\\overline{y}$ [$\mu$m]")
        _pl.xlabel("$S$ [m]")

        _pymad8.Plot.setCallbacks(figure,ax0,[ax1,ax2],self.mad8Twiss)
        
        _pl.savefig("mad8bdsim_mean.pdf")

    
    def plotBeta(self) : 
        figure = _plt.figure(figsize=(11.6, 7.2))
        gs  = _plt.GridSpec(3,1,height_ratios=[1,3,3])
        ax0 = figure.add_subplot(gs[0],projection="_My_Axes")
        _pymad8.Plot.drawMachineLattice(self.mad8Comm,self.mad8Twiss)
        
        ax1 = _pl.subplot(gs[1])
        ax1.set_autoscale_on(True)
        _pl.plot(self.mad8Envel.getColumn('suml'),_pl.sqrt(self.mad8Twiss.getColumn('betx')),"+-")
        _pl.plot(self.bdsimOptics['S'],_pl.sqrt(self.bdsimOptics['Beta_x']),"+--")
        _pl.ylabel("$\sqrt\\beta_x$ [m]")

        ax2 = _pl.subplot(gs[2])
        ax2.set_autoscale_on(True)
        _pl.plot(self.mad8Envel.getColumn('suml'),_pl.sqrt(self.mad8Twiss.getColumn('bety')),"+-")
        _pl.plot(self.bdsimOptics['S'],_pl.sqrt(self.bdsimOptics['Beta_y']),"+--")
        _pl.ylabel("$\sqrt\\beta_y$ [m]")
        _pl.xlabel("$S$ [m]")

        _pymad8.Plot.setCallbacks(figure,ax0,[ax1,ax2],self.mad8Twiss)
        
        _pl.savefig("mad8bdsim_beta.pdf")
    

    def plotSurvey(self, mad8SurveyFileName, bdsimSurveyFileName) :
        # load bdsim survey
        fs = _pybdsim.Data.Load(bdsimSurveyFileName)

        # load mad8 survey
        rs = _pymad8.Mad8.OutputReader()    
        [common, mad8Survey] = rs.readFile(mad8SurveyFileName,"survey")
        
        figure = _plt.figure(figsize=(11.6, 7.2))
        gs  = _plt.GridSpec(3,1,height_ratios=[1,3,3])
        ax0 = figure.add_subplot(gs[0],projection="_My_Axes")
        _pymad8.Plot.drawMachineLattice(self.mad8Comm,self.mad8Twiss)

        ax1 = _pl.subplot(gs[1])
        _pl.plot(mad8Survey.getColumn('suml'), mad8Survey.getColumn('x'),"+-", label = "MAD8")
        _pl.plot(fs.SStart(),fs.X(),"+--", label = "BDSIM")
        #_pl.xlim(0,max(mad8Survey.getColumn('suml')))
        _pl.ylabel("$X$ [m]")

        ax2 = _pl.subplot(gs[2])
        _pl.plot(mad8Survey.getColumn('suml'),mad8Survey.getColumn('y'),"+-", label = "MAD8")
        _pl.plot(fs.SStart(),fs.Y(),"+--", label = "BDSIM")
        #_pl.xlim(0,max(mad8Survey.getColumn('suml')))
        _pl.ylabel("$Y$ [m]")
        _pl.xlabel("$S$ [m]")

        _pl.legend(loc=0)
        _pl.subplots_adjust(hspace=0.25,top=0.94,left=0.1,right=0.92)

        _pymad8.Plot.setCallbacks(figure,ax0,[ax1,ax2],self.mad8Twiss)
        
        _pl.savefig("mad8bdsim_survey.pdf")


    def plotDispersion(self) :
        figure = _plt.figure(figsize=(11.6, 7.2))
        gs  = _plt.GridSpec(3,1,height_ratios=[1,3,3])
        ax0 = figure.add_subplot(gs[0],projection="_My_Axes")
        _pymad8.Plot.drawMachineLattice(self.mad8Comm,self.mad8Twiss)

        ax1 = _pl.subplot(gs[1])
        ax1.set_autoscale_on(True)
        _pl.plot(self.mad8Envel.getColumn('suml'),self.mad8Twiss.getColumn('dx'),"+-",label="MAD8")
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Disp_x'],"+--",label="BDSIM") # 4000 1/250*1e6
        _pl.ylabel("$\\eta_x$ [m]")

        ax2 = _pl.subplot(gs[2])
        ax2.set_autoscale_on(True)
        _pl.plot(self.mad8Envel.getColumn('suml'),self.mad8Twiss.getColumn('dy'),"+-",label="MAD8")
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Disp_y'],"+--",label="BDSIM") # 4000 1/250*1e6
        _pl.ylabel("$\\eta_y$ [m]")
        _pl.xlabel("$S$ [m]")
        
        _pl.legend(loc=0)
        
        _pymad8.Plot.setCallbacks(figure,ax0,[ax1,ax2],self.mad8Twiss)

        _pl.savefig("mad8bdsim_eta.pdf")
        

    def plotDispersionPrim(self) :
        figure = _plt.figure(figsize=(11.6, 7.2))
        gs  = _plt.GridSpec(3,1,height_ratios=[1,3,3])
        ax0 = figure.add_subplot(gs[0],projection="_My_Axes")
        _pymad8.Plot.drawMachineLattice(self.mad8Comm,self.mad8Twiss)

        ax1 = _pl.subplot(gs[1])
        ax1.set_autoscale_on(True)
        _pl.plot(self.mad8Envel.getColumn('suml'),self.mad8Twiss.getColumn('dpx'),"+-",label="MAD8")
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Disp_xp'],"+--",label="BDSIM") # 4000 1/250*1e6
        _pl.ylabel("$\\eta^{'}_{x}$ [rad]")

        ax2 = _pl.subplot(gs[2])
        ax2.set_autoscale_on(True)
        _pl.plot(self.mad8Envel.getColumn('suml'),self.mad8Twiss.getColumn('dpy'),"+-",label="MAD8")
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Disp_yp'],"+--",label="BDSIM") # 4000 1/250*1e6
        _pl.ylabel("$\\eta^{'}_{y}$ [rad]")
        _pl.xlabel("$S$ [m]")
        
        _pl.legend(loc=0)
        
        _pymad8.Plot.setCallbacks(figure,ax0,[ax1,ax2],self.mad8Twiss)

        _pl.savefig("mad8bdsim_etaprim.pdf")

        
    def plotEmittance(self) :
        figure = _plt.figure(figsize=(11.6, 7.2))
        gs  = _plt.GridSpec(3,1,height_ratios=[1,3,3])
        ax0 = figure.add_subplot(gs[0],projection="_My_Axes")
        _pymad8.Plot.drawMachineLattice(self.mad8Comm,self.mad8Twiss)

        ax1 = _pl.subplot(gs[1])
        ax1.set_autoscale_on(True)
        #_pl.plot(self.mad8Envel.getColumn('suml'), _np.zeros(len(self.mad8Envel.getColumn('suml'))),"+-",label="MAD8")
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Emitt_x'],"+--",label="BDSIM") # 4000 1/250*1e6
        _pl.ylabel("$\\epsilon_{x}$ [m]")

        ax2 = _pl.subplot(gs[2])
        ax2.set_autoscale_on(True)
        #_pl.plot(self.mad8Envel.getColumn('suml'), _np.zeros(len(self.mad8Envel.getColumn('suml'))),"+-",label="MAD8")
        _pl.plot(self.bdsimOptics['S'],self.bdsimOptics['Emitt_y'],"+--",label="BDSIM") # 4000 1/250*1e6
        _pl.ylabel("$\\epsilon_{y}$ [m]")
        _pl.xlabel("$S$ [m]")
        
        _pl.legend(loc=0)
        
        _pymad8.Plot.setCallbacks(figure,ax0,[ax1,ax2],self.mad8Twiss)

        _pl.savefig("mad8bdsim_emitt.pdf")


