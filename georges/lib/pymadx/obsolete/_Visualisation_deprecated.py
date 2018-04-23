
#the Tfs data class
from Tfs import Tfs

import pymadx as _pymadx
import pylab as _pylab
import numpy as _np
import matplotlib.patches as _patches
import matplotlib.pyplot   as _plt


class Visualise:
    def __init__(self,tfsobject):
        if type(tfsobject) == str:
            self.tfsfile = Tfs(tfsobject)
        elif type(tfsobject) == Tfs:
            self.tfsfile = Tfs(tfsobject)
        else:
            raise IOError("Visualise> Unkown input type "+str(tfsobject))

    def _ExtractRelevantData(self):
        self.data = {}
        self.data['S'] = _np.array(self.tfs.ColumnByIndex('S'))
        self._kw   = _np.array(self.tfs.ColumnByIndex('KEYWORD'))
        self._l    = _np.array(self.tfs.ColumnByIndex('L'))
        self._betx = _np.array(self.tfs.ColumnByIndex('BETX'))
        self._bety = _np.array(self.tfs.ColumnByIndex('BETY'))
        self._alpx = _np.array(self.tfs.ColumnByIndex('ALPX'))
        self._alpy = _np.array(self.tfs.ColumnByIndex('ALPX'))

class Visualisation(object) : 
    def __init__(self, tfsFile = "../madx/ff_twiss", surveyFile = "") : 
        self.tfsFile = tfsFile
        self.tfs  = _pymadx.Tfs(tfsFile) 
        self.s    = _np.array(self.tfs.ColumnByIndex('S'))
        self.kw   = _np.array(self.tfs.ColumnByIndex('KEYWORD'))
        self.l    = _np.array(self.tfs.ColumnByIndex('L'))
        self.betx = _np.array(self.tfs.ColumnByIndex('BETX'))
        self.bety = _np.array(self.tfs.ColumnByIndex('BETY'))
#        self.alpx = _np.array(self.tfs.ColumnByIndex('ALPX'))
#        self.alpy = _np.array(self.tfs.ColumnByIndex('ALPX'))

# mangets inc multipoles         
        self.k0l  = _np.array(self.tfs.ColumnByIndex('K0L'))
        self.k1l  = _np.array(self.tfs.ColumnByIndex('K1L'))
        self.k2l  = _np.array(self.tfs.ColumnByIndex('K2L'))
        self.k3l  = _np.array(self.tfs.ColumnByIndex('K3L'))
        self.k4l  = _np.array(self.tfs.ColumnByIndex('K4L'))
        self.k5l  = _np.array(self.tfs.ColumnByIndex('K5L'))

                
    def draw(self) : 
        fig = _pylab.figure(1) 
        gs = _pylab.GridSpec(2,1,height_ratios=(1,4));
        
        axMachine = _pylab.subplot(gs[0])
        axMachine.axes.get_yaxis().set_visible(False)
        self.drawMachine()

        axOptics = _pylab.subplot(gs[1])
        self.plotOptics()
        
        def machineXlim(ax) : 
            axMachine.set_autoscale_on(False)
            axOptics.set_xlim(axMachine.get_xlim())

        def click(a) : 
            if a.button == 3 : 
                print self.tfs.NameFromNearestS(a.xdata)

        axMachine.callbacks.connect('xlim_changed', machineXlim) 
        fig.canvas.mpl_connect('button_press_event', click) 
            
    def drawMachine(self) : 
        ############## draw bend
        def drawBend(s,l) : 
            br = _patches.Rectangle((s,-0.1),l,0.2,color='g')
            ax = _plt.gca()
            ax.add_patch(br)
        ############## draw quad
        def drawQuad(s,l,k) : 
            if k > 0 :
                qr = _patches.Rectangle((s,0),l,0.2,color='r')
            elif k < 0 : 
                qr = _patches.Rectangle((s,-0.2),l,0.2,color='b')
            else : 
                qr = _patches.Rectangle((s,-0.1),l,0.2,color='g')                

            ax = _plt.gca()
            ax.add_patch(qr)
        ############## draw sextupole
        def drawSext(s,l,k) : # 6-a-gon
            pass 
        ############## draw octupole 
        def drawOct(s,l,k) : # Octagon
            pass 
        ############## draw multipole    
        def drawMulti(s,l) : # ellipse
            pass
        def drawMarker(s) :         
            _plt.plot([s,s],[0.4,0.6],"-")
#            _pylab.axvline(s)
            

        ############## draw monitor
        def drawMonitor(s) : 
            _plt.plot([s,s],[0.65,0.85],"-")            
            pass
        ############## draw collimator
        def drawCollimator(s,l) : 
            br = _patches.Rectangle((s,-0.6),l,0.2,color='k')
            ax = _plt.gca()
            ax.add_patch(br)
            
        # plot beam line 
        z = _pylab.zeros(self.s.shape)
        _pylab.plot(self.s,z,'k--')
        _pylab.ylim(-1.0,1.0)

        # plot lines for monitors
        _pylab.plot([0,max(self.s)],[0.4,0.4],'k')
        _pylab.plot([0,max(self.s)],[0.6,0.6],'k')

        # plot lines for markers 
        _pylab.plot([0,max(self.s)],[0.65,0.65],'k')
        _pylab.plot([0,max(self.s)],[0.85,0.85],'k')


        # plot lines for collimators
        _pylab.plot([0,max(self.s)],[-0.4,-0.4],'k')        
        _pylab.plot([0,max(self.s)],[-0.6,-0.6],'k')  

    def DrawMachineSimple(self) : 
        def DrawBend(s,l) : 
            br = _patches.Rectangle((s,-0.1),l,0.2,color='g')
            ax = _plt.gca()
            ax.add_patch(br)
        def DrawQuad(s,l,k) : 
            if k > 0 :
                qr = _patches.Rectangle((s,0),l,0.2,color='r')
            elif k < 0 : 
                qr = _patches.Rectangle((s,-0.2),l,0.2,color='b')
            else : 
                qr = _patches.Rectangle((s,-0.1),l,0.2,color='#B2B2B2') #a nice grey in hex
            ax = _plt.gca()
            ax.add_patch(qr)
        def DrawSext(s,l,k) : # 6-a-gon
            pass 
        ############## Draw octupole 
        def DrawOct(s,l,k) : # Octagon
            pass 
        ############## Draw multipole    
        def DrawMulti(s,l) : # ellipse
            pass
        #def DrawMarker(s) :         
        #    _plt.plot([s,s],[0.4,0.6],"-")
#            _pylab.axvline(s)
            

        ############## Draw monitor
        #def DrawMonitor(s) : 
        #    _plt.plot([s,s],[0.65,0.85],"-")            
        #    pass
        ############## Draw collimator
        def DrawCollimator(s,l) : 
            br = _patches.Rectangle((s,-0.2),l,0.4,color='k')
            ax = _plt.gca()
            ax.add_patch(br)
            
        # plot beam line 
        z = _pylab.zeros(self.s.shape)
        _pylab.plot(self.s,z,'k-')
        _pylab.ylim(-1.0,1.0)

        # plot lines for monitors
        #_pylab.plot([0,max(self.s)],[0.4,0.4],'k')
        #_pylab.plot([0,max(self.s)],[0.6,0.6],'k')

        # plot lines for markers 
        #_pylab.plot([0,max(self.s)],[0.65,0.65],'k')
        #_pylab.plot([0,max(self.s)],[0.85,0.85],'k')


        # plot lines for collimators
        #_pylab.plot([0,max(self.s)],[-0.4,-0.4],'k')        
        #_pylab.plot([0,max(self.s)],[-0.6,-0.6],'k')  


        # loop over elements and Draw on beamline         
        for i in range(0,len(self.s),1) : 

            if self.kw[i] == 'QUADRUPOLE' : 
                DrawQuad(self.s[i], self.l[i], self.k1l[i])
            elif self.kw[i] == 'RBEND' : 
                DrawBend(self.s[i], self.l[i])
            elif self.kw[i] == 'SBEND' : 
                DrawBend(self.s[i], self.l[i])
            elif self.kw[i] == 'RCOLLIMATOR' : 
                DrawCollimator(self.s[i],self.l[i])
            elif self.kw[i] == 'ECOLLIMATOR' : 
                DrawCollimator(self.s[i],self.l[i])


        
    def plotOptics(self) : 
        _pylab.plot(self.s,_np.sqrt(self.betx),'b-+')#,label="$\\sqrt{\beta_x}$")
        _pylab.plot(self.s,_np.sqrt(self.bety),'g-+')#,label="$\\sqrt{\beta_y}$")
        #_pylab.xlabel("$S$~(m)")
        #_pylab.ylabel("$\sqrt{\beta_{x,y}}$~(${\rm{m}}^2$)")
        

        _pylab.legend(loc=2)

    def plotBeamSize(self, nemitx = 0.0, nemity = 0.0) :
        pass

