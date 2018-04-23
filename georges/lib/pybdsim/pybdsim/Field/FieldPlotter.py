import matplotlib.pyplot as _plt
import numpy as _np

#import pybdsim

class FourDData(object):
    def __init__(self, filename, xind=0, yind=1, zind=2, tind=3):
        d = pybdsim.Field.Load(filename)
        
        # '...' fills in unknown number of dimensions with ':' meaning
        # all of that dimension
        if (xind >= 0):
            self.x  = d[..., xind].flatten()
        if (yind >= 0):
            self.y  = d[..., yind].flatten()
        if (zind >= 0):
            self.z  = d[..., zind].flatten()
        if (tind >= 0):
            self.t  = d[..., tind].flatten()

        # index from end as we don't know the dimensionality
        self.fx = d[...,-3].flatten()
        self.fy = d[...,-2].flatten()
        self.fz = d[...,-1].flatten()

        self.mag = _np.sqrt(self.fx**2 + self.fy**2 + self.fz**2)

class ThreeDData(FourDData):
    def __init__(self, filename):
        FourDData.__init__(self, filename, tind=-1)

class TwoDData(FourDData):
    def __init__(self, filename):
        FourDData.__init__(self, filename, tind=-1, zind=-1)

class OneDData(FourDData):
    def __init__(self, filename):
        FourDData.__init__(self, filename, tind=-1, zind=-1, yind=-1)

def _Niceties():
    _plt.xlabel('X (cm)')
    _plt.ylabel('Y (cm)')
    _plt.colorbar()
    _plt.tight_layout()
    _plt.axes().set_aspect('equal', 'datalim')

def Plot2DXY(filename, scale=None):
    d = TwoDData(filename)
    _plt.figure()
    _plt.quiver(d.x,d.y,d.fx,d.fy,d.mag,cmap=_plt.cm.magma,pivot='mid',scale=scale)
    _Niceties()

def Plot3DXY(filename, scale=None):
    d = ThreeDData(filename)
    _plt.figure()
    _plt.quiver(d.x,d.y,d.fx,d.fy,d.mag,cmap=_plt.cm.magma,pivot='mid',scale=scale)
    _Niceties()

def Plot3DXZ(filename, scale=None):
    d = ThreeDData(filename)
    _plt.figure()
    _plt.quiver(d.x,d.z,d.fx,d.fz,d.mag,cmap=_plt.cm.magma,pivot='mid',scale=scale)
    _Niceties()
