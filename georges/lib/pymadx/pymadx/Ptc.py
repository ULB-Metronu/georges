"""
Classes to handle PTC runs and data.
"""

import numpy as _np
import re as _re
try:
    import matplotlib.pyplot as _plt
except ImportError:
    pass
from numpy.random import multivariate_normal as _multivariate_normal

class Inray(object):
    """
    Class for a madx ptc input ray
    x  : horizontal position [m]
    px : horizontal canonical momentum p_x/p_0 
    y  : vertical position [m]
    py : vertical canonical momentum p_y/p_0
    t  : c*(t-t0) [m]
    pt : (delta-E)/(pc)

    use str(Inray) to get the representation for file writing
    """
    def __init__(self, x=0.0, px=0.0, y=0.0, py=0.0, t=0.0, pt=0.0):
        self.x  = x
        self.px = px
        self.y  = y
        self.py = py
        self.t  = t
        self.pt = pt

    def __repr__(self):
        s =  'ptc_start'
        s += ', x='  + str(self.x)
        s += ', px=' + str(self.px)
        s += ', y='  + str(self.y)
        s += ', py=' + str(self.py)
        s += ', t='  + str(self.t)
        s += ', pt=' + str(self.pt)
        s += ';\n'
        return s

class Inrays(list):
    """
    Class based on python list for Inray class
    """
    def __init__(self):
        list.__init__(self)
        variables = ['X','PX','Y','PY','T','PT']
        for v in variables:
            self._AddMethod(v)

    def AddParticle(self,x=0.0,px=0.0,y=0.0,py=0.0,t=0.0,pt=0.0):
        self.append(Inray(x,px,y,py,t,pt))

    def Clear(self):
        del self[:]

    def Write(self,filename):
        WriteInrays(filename,self)

    def Plot(self):
        PlotInrays(self)

    def _AddMethod(self, variablename):
        """This is used to easily and dynamically add a getter function for a variable name."""
        def GetAttribute():
            return _np.array([getattr(p,str.lower(variablename)) for p in self])
        setattr(self,variablename,GetAttribute)

    def Statistics(self):
        print('TBC - will return various moments')
        
def LoadInrays(fileName): 
    """Load input rays from file
    fileName : inrays.madx 
    return   : Inrays instance""" 
    i = Inrays()
    
    # open file 
    f = open(fileName) 
    for l in f : 
        inre_x  = _re.search('\s*x\s*=\s*([0-9.eE+-]+)\s*',l)
        inre_px = _re.search('\s*px\s*=\s*([0-9.eE+-]+)\s*',l)
        inre_y  = _re.search('\s*y\s*=\s*([0-9.eE+-]+)\s*',l)
        inre_py = _re.search('\s*py\s*=\s*([0-9.eE+-]+)\s*',l)
        inre_t  = _re.search(' t\s*=\s*([0-9.eE+-]+)\s*',l)
        inre_pt = _re.search('\s*pt\s*=\s*([0-9.eE+-]+)\s*',l)

        if inre_x : 
            x = float(inre_x.group(1))
        else :
            x = 0.0
        if inre_px : 
            px = float(inre_px.group(1))
        else : 
            px = 0.0 
        if inre_y : 
            y  = float(inre_y.group(1))
        else : 
            y  = 0.0
        if inre_py :             
            py = float(inre_py.group(1))
        else : 
            py = 0.0
        if inre_t :             
            t  = float(inre_t.group(1))
        else : 
            t = 0.0

        pt = float(inre_pt.group(1))

        i.AddParticle(x,px,y,py,t,pt)

    print('LoadInrays> Loaded ',len(i))
    return i
  
def WriteInrays(fileName, inrays):
    print('pymadx.Ptc> WriteInrays - inrays written to: ',fileName)
    f = open(fileName, 'w') 
    for particle in inrays:
        f.write(str(particle))
    f.close()
   
def PlotInrays(i): 
    """Plot Inrays instance, if input is a sting the instance is created from the file"""    

    if type(i) == str : 
        i = LoadInrays(i)
        
    f = _plt.figure(1) 
    f.clf()
    
    ax = f.add_subplot(3,2,1)    
    ax.hist(i.X(),50,histtype='step',label='x')
    ax.text(0.05,0.85,'X',transform=ax.transAxes)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax = f.add_subplot(3,2,2)    
    ax.hist(i.PX(),50,histtype='step',label='px')
    ax.text(0.05,0.85,'PX',transform=ax.transAxes)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    ax = f.add_subplot(3,2,3)    
    ax.hist(i.Y(),50,histtype='step',label='y')
    ax.text(0.05,0.85,'Y',transform=ax.transAxes)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax = f.add_subplot(3,2,4)    
    ax.hist(i.PY(),50,histtype='step',label='py')
    ax.text(0.05,0.85,'PY',transform=ax.transAxes)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    ax = f.add_subplot(3,2,5)    
    ax.hist(i.T(),50,histtype='step',label='t')
    ax.text(0.05,0.85,'T',transform=ax.transAxes)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax = f.add_subplot(3,2,6)    
    ax.hist(i.PT(),50,histtype='step',label='pt')
    ax.text(0.05,0.85,'PT',transform=ax.transAxes)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    _plt.subplots_adjust(hspace=0.35,wspace=0.15,top=0.98,right=0.98,left=0.05)
  
class GaussGenerator(object): 
    """Simple ptx inray file generator"""
    def __init__(self,
                 gemx = 1e-10, betax = 0.1, alfx = 0.0 , 
                 gemy = 1e-10, betay = 0.1, alfy = 0.0,
                 sigmat = 1e-12, sigmapt= 1e-12): 
        """Simple gaussian beam
        gemx   : x geometric emittance 
        betax  : x beta function
        alfx   : x alpha function
        gemy   : y geometric emittance
        betay  : y beta function
        alfy   : y alpha function 
        sigmat : gaussian spread in time around the reference time 
        sigmapt: gaussian spread in relative energy 
        """
        
        self.gemx    = gemx
        self.betax   = betax
        self.alfx    = alfx

        self.gemy    = gemy
        self.betay   = betay
        self.alfy    = alfy 

        self.sigmat  = sigmat 
        self.sigmapt = sigmapt

        self.gamx    = (1.0+self.alfx**2)/self.betax    
        self.gamy    = (1.0+self.alfy**2)/self.betay

        # Generate sigma matrix
        self.means  = _np.zeros(6)
        self.sigmas = _np.zeros((6,6)) 
        
        self.sigmas[0][0] =  self.gemx*self.betax
        self.sigmas[0][1] = -self.gemx*self.alfx
        self.sigmas[1][0] = -self.gemx*self.alfx
        self.sigmas[1][1] =  self.gemx*self.gamx
        self.sigmas[2][2] =  self.gemy*self.betay
        self.sigmas[2][3] = -self.gemy*self.alfy
        self.sigmas[3][2] = -self.gemy*self.alfy
        self.sigmas[3][3] =  self.gemy*self.gamy
        self.sigmas[4][4] =  self.sigmat**2
        self.sigmas[5][5] =  self.sigmapt**2
        
    def __repr__(self) : 
        s = 'ex : '+str(self.gemx)+' bx : '+str(self.betax)+' ax : '+str(self.alfx)+' gx : '+str(self.gamx)+'\n'
        s+= 'ey : '+str(self.gemy)+' bx : '+str(self.betay)+' ax : '+str(self.alfy)+' gy : '+str(self.gamy)+'\n'
        s+= 'sT : '+str(self.sigmat)+' spt : '+str(self.sigmapt)
        return s

    def Generate(self, nToGenerate=1000, fileName='inrays.madx'): 
        """ returns an Inrays structure""" 
        i = Inrays()
        
        for c in range(0,nToGenerate,1):
            rv = _multivariate_normal(self.means,self.sigmas,1)[0]
            i.AddParticle(rv[0],rv[1],rv[2],rv[3],rv[4],rv[5]) 

        WriteInrays(fileName,i)
        #return i

class FlatGenerator(object):
    """Simple ptc inray file generator - even distribution"""
    def __init__(self,
                 mux =0.0, widthx =1e-3,
                 mupx=0.0, widthpx=1e-3,
                 muy =0.0, widthy =1e-3,
                 mupy=0.0, widthpy=1e-3):
        self.mux     = mux
        self.muy     = muy
        self.widthx  = widthx
        self.widthy  = widthy
        self.mupx    = mupx
        self.mupy    = mupy
        self.widthpx = widthpx
        self.widthpy = widthpy

    def __repr__(self):
        s = '' #TBC
        return s

    def Generate(self,nToGenerate=100, fileName='inrays.madx'):
        """ returns an Inrays structure""" 
        i = Inrays()

        nd = 0.0
        if self.widthx > 0:
            nd +=1
        if self.widthy > 0:
            nd +=1
        if self.widthpx > 0:
            nd +=1
        if self.widthpy > 0:
            nd +=1
        nperdim = _np.ceil(nToGenerate**(1/nd))
        print("FlatGenerator> making array square - there'll be ",nperdim**nd,'particles')

        xmin = self.mux - self.widthx/2.0
        xmax = self.mux + self.widthx/2.0
        ymin = self.muy - self.widthy/2.0
        ymax = self.muy + self.widthy/2.0
        pxmin = self.mupx - self.widthpx/2.0
        pxmax = self.mupx + self.widthpx/2.0
        pymin = self.mupy - self.widthpy/2.0
        pymax = self.mupy + self.widthpy/2.0

        def GetX():
            return _np.linspace(ymin,ymax,nperdim)

        def GetY():
            return _np.linspace(xmin,xmax,nperdim)

        def GetXp():
            return _np.linspace(pymin,pymax,nperdim)

        def GetYp():
            return _np.linspace(pxmin,pxmax,nperdim)

        #we don't actually use 'd' here but just take advantage of list comprehensions
        d = [i.AddParticle(x=x,px=px,y=y,py=py) for x in GetX() for px in GetXp() for y in GetY() for py in GetYp()]

        WriteInrays(fileName,i)
        #return i
