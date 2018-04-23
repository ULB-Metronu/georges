import pymadx.Ptc
import pymadx.Beam
import pymadx.Builder
import pymadx.Tfs
import pybdsim.Beam
import pybdsim.Builder
import pybdsim.Data
import os as _os
import matplotlib.pyplot as _plt
import robdsim
import numpy as _np

class Test:
    def __init__(self,type_,foldername=None,particle="e-",energy=1.0,distribution='flat',nparticles=10,length=1.0,**kwargs): 
        """
        Tracking test class

        type = 'drift' | 'quadrupole' | 'sextupole' | 'sbend' | 'solenoid'
        distribution = 'flat' | 'gaussian'
        """
        self.type_        = type_
        self.filename     = self.type_
        self.foldername   = foldername
        self.ptcfilename  = 'inrays.madx'
        if self.foldername != None:
            self.usingfolder = True
            self.filepath = self.foldername+'/'+self.filename
            self.ptcfilepath =  self.foldername+'/'+self.ptcfilename
            _os.system("mkdir -p " + self.foldername)
        else:
            self.usingfolder = False
            self.filepath = self.filename
            self.ptcfilepath = self.ptcfilename
        self.particle     = particle
        self.energy       = energy
        self.distribution = distribution
        self.distrkwargs  = {}
        self.nparticles   = nparticles
        self.length       = length
        self.kwargs       = kwargs
        self.figureNr     = 1729 # arbitrary number where figure start from
        
    def CleanMakeRun(self):
        self.Clean()
        self.Make()
        self.Execute()
        self.Compare()
    
    def Clean(self):        
        _os.system("rm -rf "+self.filepath+"*")
        _os.system("rm -rf "+self.foldername+"/output*")
        _os.system("rm -rf "+self.foldername+"/Maxwellian_bend_for_ptc.txt trackone inrays.madx")
        _os.system("rm -rf "+self.foldername+"/test*")
        _os.system("rm -rf "+self.foldername+"/*.log")
        _os.system("rm -rf "+self.foldername+"/*.dat")        
        _os.system("rm -rf "+self.foldername+"/trackone")
        _os.system("rm -rf "+self.foldername+"/inrays.madx")

        # clean and close figures (8 figures in total)
        for i in range(9):
            _plt.close(self.figureNr+i)
        
    def ChangeDistribution(self,distribution='flat',nparticles=10,**kwargs):
        """
        'flat'
        kwargs: mux=0.0, widthx=1e-05, mupx=0.0, widthpx=1e-05, muy=0.0, 
                widthy=0.001, mupy=0.0, widthpy=0.001
        'gauss'
        kwargs: gemx=1e-10, betax=0.1, alfx=0.0, gemy=1e-10, betay=0.1, 
                alfy=0.0, sigmat=1e-12, sigmapt=1e-12
        """
        self.distribution = distribution
        self.distrkwargs  = kwargs
        self.nparticles   = nparticles
        
    def Make(self):
        #type_='drift', foldername=None, particle='e-', energy=1.0,**kwargs) : 
        print 'Test> Element type:         ',self.type_
        print 'Test> Destination filepath: ',self.filepath
        print 'Test> kwargs: ',
        for k in self.kwargs : 
            print k+'='+str(self.kwargs[k]),

        if self.distribution == 'flat':
            ptc = pymadx.Ptc.FlatGenerator(**self.distrkwargs)
        elif self.distribution == 'gauss':
            ptc = pymadx.Ptc.GaussGenerator(**self.distrkwargs)
        
        ptc.Generate(self.nparticles,self.ptcfilepath)
        
        bb  = pybdsim.Beam.Beam(self.particle,self.energy,'ptc')
        xb  = pymadx.Beam(self.particle,self.energy,'ptc')
        
        bb.SetDistribFileName(self.ptcfilename) 
        xb.SetDistribFileName(self.ptcfilename) 
        
        bm  = pybdsim.Builder.Machine()
        xm  = pymadx.Builder.Machine()
        
        bm.AddBeam(bb)
        xm.AddBeam(xb)
    
        if self.type_ == 'drift' :
            name = 'd1'
            bm.AddDrift(name,length=self.length,**self.kwargs)
            xm.AddDrift(name,length=self.length,**self.kwargs)
        elif self.type_ == 'sbend':
            name = 'sb1'
            bm.AddDipole(name,length=self.length,**self.kwargs)
            xm.AddDipole(name,length=self.length,**self.kwargs)
        elif self.type_ == 'rbend':
            name = 'rb1'
            bm.AddDipole(name,category='rbend',length=self.length,**self.kwargs)
            xm.AddDipole(name,category='rbend',length=self.length,**self.kwargs)
        elif self.type_ == 'hkicker':
            name = 'hk1'
            bm.AddHKicker(name, length=self.length, **self.kwargs)
            xm.AddHKicker(name, length=self.length, **self.kwargs)
        elif self.type_ == 'vkicker':
            name = 'vk1'
            bm.AddVKicker(name, length=self.length, **self.kwargs)
            xm.AddVKicker(name, length=self.length, **self.kwargs)
        elif self.type_ == 'quadrupole' :
            name = 'q1'
            bm.AddQuadrupole(name,length=self.length,**self.kwargs)
            xm.AddQuadrupole(name,length=self.length,**self.kwargs)
        elif self.type_ == 'sextupole' :
            name = 's1'
            bm.AddSextupole(name,length=self.length,**self.kwargs)
            xm.AddSextupole(name,length=self.length,**self.kwargs)
        elif self.type_ == 'solenoid':
            name = 'sl1'
            bm.AddSolenoid(name,length=self.length,**self.kwargs)
            xm.AddSolenoid(name,length=self.length,**self.kwargs)

        bm.AddMarker("theend") # Need a post element marker to sample at, only for bdsim
        
        bm.AddSampler('theend')
        xm.AddSampler('all')

        bm.Write(self.filepath)
        xm.Write(self.filepath)
            
    def Execute(self):
        if self.usingfolder:
            _os.chdir(self.foldername)
        
        _os.system("madx < "+self.filename+".madx > madx.log")
        _os.system("bdsim --file="+self.filename+".gmad --batch --output=combined --outfile='test' > bdsim.log")
        
        if self.usingfolder:
            _os.chdir("../")

    def Compare(self, addPrimaries=True):

        if self.usingfolder:
            _os.chdir(self.foldername)

        bdsimprim = pybdsim.Data.Load("test/test.primaries.txt")
        Bx0 = bdsimprim.X()
        By0 = bdsimprim.Y()
        Bxp0 = bdsimprim.Xp()
        Byp0 = bdsimprim.Yp()
        self.bdsimprimaries = {'x':Bx0,'y':By0,'xp':Bxp0,'yp':Byp0}
        
        bdsim = pybdsim.Data.Load("test/test.txt")
        Bx = bdsim.X()
        By = bdsim.Y()
        Bxp = bdsim.Xp()
        Byp = bdsim.Yp()
        self.bdsimoutput = {'x':Bx,'y':By,'xp':Bxp,'yp':Byp}

        madx = pymadx.Tfs("trackone")
        madx = madx.GetSegment(madx.nsegments) #get the last 'segment' / sampler
        Mx = madx.GetColumn('X')*1e6 # convert from m to um
        My = madx.GetColumn('Y')*1e6 
        Mxp = madx.GetColumn('PX')
        Myp = madx.GetColumn('PY')
        self.ptcoutput = {'x':Mx,'y':My,'xp':Mxp,'yp':Myp}

        fresx  = _np.nan_to_num(Mx - Bx)
        fresy  = _np.nan_to_num(My - By)
        fresx  = _np.nan_to_num(fresx / Mx) #protect against nans for 0 diff
        fresy  = _np.nan_to_num(fresy / My)
        fresxp = _np.nan_to_num(Mxp - Bxp)
        fresyp = _np.nan_to_num(Myp - Byp)
        fresxp = _np.nan_to_num(fresxp / Mxp)
        fresyp = _np.nan_to_num(fresyp / Myp)
        self.residuals = {'x':fresx,'y':fresy,'xp':fresxp,'yp':fresyp}
        
        # 2d plots
        #X vs Y
        _plt.figure(self.figureNr)
        _plt.clf()
        _plt.plot(Mx,My,'b.',label='PTC')
        _plt.plot(Bx,By,'g.',label='BDSIM')
        if addPrimaries:
            _plt.plot(Bx0,By0,'r.',label='BDSIM prim')
        _plt.legend()
        _plt.xlabel(r"x ($\mu$m)")
        _plt.ylabel(r"y ($\mu$m)")
        _plt.title(self.type_)
        _plt.savefig(self.type_+'_xy.pdf')
        _plt.savefig(self.type_+'_xy.png')

        #XP vs YP
        _plt.figure(self.figureNr+1)
        _plt.clf()
        _plt.plot(Mxp,Myp,'b.',label='PTC')
        _plt.plot(Bxp,Byp,'g.',label='BDSIM')
        if addPrimaries:
            _plt.plot(Bxp0,Byp0,'r.',label='BDSIM prim')
        _plt.legend()
        _plt.xlabel(r"x' ($\mu$m)")
        _plt.ylabel(r"y' ($\mu$m)")
        _plt.title(self.type_)
        _plt.savefig(self.type_+'_xpyp.pdf')
        _plt.savefig(self.type_+'_xpyp.png')

        #X vs XP
        _plt.figure(self.figureNr+2)
        _plt.clf()
        _plt.plot(Mx,Mxp,'b.',label='PTC')
        _plt.plot(Bx,Bxp,'g.',label='BDSIM')
        if addPrimaries:
            _plt.plot(Bx0,Bxp0,'r.',label='BDSIM prim')
        _plt.legend()
        _plt.xlabel(r"x ($\mu$m)")
        _plt.ylabel(r"x' (rad)")
        _plt.title(self.type_)
        _plt.savefig(self.type_+'_xxp.pdf')
        _plt.savefig(self.type_+'_xxp.png')

        #Y vs YP
        _plt.figure(self.figureNr+3)
        _plt.clf()
        _plt.plot(My,Myp,'b.',label='PTC')
        _plt.plot(By,Byp,'g.',label='BDSIM')
        if addPrimaries:
            _plt.plot(By0,Byp,'r.',label='BDSIM prim')
        _plt.legend()
        _plt.xlabel(r"y ($\mu$m)")
        _plt.ylabel(r"y' (rad)")
        _plt.title(self.type_)
        _plt.savefig(self.type_+'_yyp.pdf')
        _plt.savefig(self.type_+'_yyp.png')

        # 1d plots
        # x comparison
        f = _plt.figure(self.figureNr+4)
        f.suptitle(self.type_)
        _plt.clf()

        ax1 = f.add_subplot(221)
        ax1.hist(Mx,color='b',label='PTC',histtype='step')
        ax1.hist(Bx,color='g',label='BDSIM',histtype='step')
        if addPrimaries:
            ax1.hist(Bx0,color='r',label='BDSIM prim',histtype='step')
        ax1.legend(fontsize='x-small',loc=0)
        ax1.set_xlabel(r"x ($\mu$m)")
        
        # y comparison
        ax2 = f.add_subplot(222)
        ax2.hist(My,color='b',label='PTC',histtype='step')
        ax2.hist(By,color='g',label='BDSIM',histtype='step')
        if addPrimaries:
            ax2.hist(By0,color='r',label='BDSIM prim',histtype='step')
        ax2.legend(fontsize='x-small',loc=0)
        ax2.set_xlabel(r"y ($\mu$m)")

        # xp comparison
        ax3 = f.add_subplot(223)
        ax3.hist(Mxp,color='b',label='PTC',histtype='step')
        ax3.hist(Bxp,color='g',label='BDSIM',histtype='step')
        if addPrimaries:
            ax3.hist(Bxp0,color='r',label='BDSIM prim',histtype='step')
        ax3.legend(fontsize='x-small',loc=0)
        ax3.set_xlabel(r"x' (rad)")

        # yp comparison
        ax4 = f.add_subplot(224)
        ax4.hist(Myp,color='b',label='PTC',histtype='step')
        ax4.hist(Byp,color='g',label='BDSIM',histtype='step')
        if addPrimaries:
            ax4.hist(Byp0,color='r',label='BDSIM prim',histtype='step')
        ax4.legend(fontsize='x-small',loc=0)
        ax4.set_xlabel(r"y' (rad)")
        
        _plt.savefig(self.type_+'_hist.pdf')
        _plt.savefig(self.type_+'_hist.png')

        # residuals in one plot
        f = _plt.figure(self.figureNr+8)
        _plt.clf()
        
        axX = f.add_subplot(221)
        axX.hist(Mx,weights=fresx,bins=100,ec='b')
        axX.set_xlabel(r'X ($\mu$m)')
        axX.set_ylabel('Fractional Residuals')
        
        axY = f.add_subplot(222)
        axY.hist(My,weights=fresy,bins=100,ec='b')
        axY.set_xlabel(r'Y ($\mu$m)')
        axY.set_ylabel('Fractional Residuals')
        
        axXp = f.add_subplot(223)
        axXp.hist(Mxp*1e3,weights=fresxp,bins=100,ec='b')
        axXp.set_xlabel('Xp (mrad)')
        axXp.set_ylabel('Fractional Residuals')

        axYp = f.add_subplot(224)
        axYp.hist(Myp*1e3,weights=fresyp,bins=100,ec='b')
        axYp.set_xlabel('Yp (mrad)')
        axYp.set_ylabel('Fractional Residuals')

        _plt.subplots_adjust(left=0.15,right=0.95,top=0.95,wspace=0.39,hspace=0.25)
        _plt.savefig(self.type_+'_residuals.pdf')
        _plt.savefig(self.type_+'_residuals.png')
        
        #emittance
        r = robdsim.robdsimOutput("test.root")
        r.CalculateOpticalFunctions("optics.dat")
        d = pybdsim.Data.Load("optics.dat")
        print 'Horizontal emittance bdsim (before,after) ',d.Emitt_x()
        print 'Vertical emittance bdsim (before,after) ',d.Emitt_y()

        if self.usingfolder:
            _os.chdir("../")
