"""
| Survey() - survey a gmad lattice, plot element coords
| Loader() - load a gmad file using the compiled bdsim parser
| GmadFile() - modify a text based gmad file
"""

import ctypes as _ctypes
import ctypes as __ctypes
from ctypes.util import find_library as _find_library
import os as _os
import numpy as _np
import matplotlib.pyplot as _plt
try:
    from StringIO import StringIO as _StringIO
except ImportError:
    from io import StringIO as _StringIO
import re as _re

class Survey(object):
    """
    Survey - load a gmad lattice and have a look

    Example:

    >>> a = Survey()
    >>> a.Load('mylattice.gmad')
    >>> a.Plot()

    """
    def __init__(self,filename=None) :
        self._z_current = 0
        self._x_current = 0
        self._theta     = _np.pi
        self._names     = []
        self._beamline  = []
        self._z_coords  = []
        self._x_coords  = []
        self._lentotal  = 0
        if filename != None:
            self.Load(filename)

    def Load(self, filename):
        print('Warning - we must check the sign of angle so dx has correct sign')
        self._file = Lattice(filename)

        for e in self._file.lattice:
            name   = e['Name']
            length = e['Length']
            angle  = e['Angle']
            self.Step(angle,length)

    def Step(self,angle,length):
        self._theta += angle
        dz    = length * _np.cos(self._theta)
        dx    = length * _np.sin(self._theta)
        self._lentotal += length

        z_new = self._z_current - dz
        x_new = self._x_current - dx
        self._beamline.append([[self._z_current,z_new],[self._x_current,x_new]])
        self._z_current = z_new
        self._x_current = x_new
        self._z_coords.append(z_new)
        self._x_coords.append(x_new)

    def FinalDiff(self):
        dz = self._z_current - self._beamline[0][0][1]
        dx = self._x_current - self._beamline[0][1][1]
        print('Final dz ',dz,'m')
        print('Final dx ',dx,'m')

    def Plot(self):
        _plt.figure()
        _plt.plot(self._z_coords,self._x_coords,'b.')
        _plt.plot(self._z_coords,self._x_coords,'b-')
        _plt.xlabel('Z (m)')
        _plt.ylabel('X (m)')
        _plt.show()

    def CompareMadX(self, fileName):
        print('NOT WORKING - TBC')
        return
        import pymadx as _pymadx

        mxs = _pymadx.Survey()
        mxs.Load(fileName)

        self.Plot()

        _plt.plot(mxs._x_coords,mxs._y_coords,'rx')
        _plt.plot(mxs._x_coords,mxs._y_coords,'r-')
        _plt.show()

    def FindClosestElement(self, coord) :
        arr = _np.arange(0,len(self._x_coords),1)
        xdiff = _np.array(self._x_coords)-coord[0]
        ydiff = _np.array(self._y_coords)-coord[1]
        d     = _np.sqrt(xdiff**2+ydiff**2)
        print(arr[d == min(d)])


def _LoadLib():
    libname = 'gmadShared'
    libpath = '/usr/local/lib/'
    parserlib = None

    test0,test1,test2,test3 = False, False, False, False

    #test0 - try general find function
    libpath0 = _find_library(libname)
    if libpath0 != None:
        parserlib = _ctypes.cdll.LoadLibrary(libpath0)
        test0 = True

    #test1 - pybdsim is being used from bdsim build dir
    _this_dir, _this_filename = _os.path.split(__file__)
    libpath1 = _this_dir+'/../../parser/'+'lib'+libname
    if _os.path.exists(libpath1+'.so'):
        parserlib = _ctypes.cdll.LoadLibrary(libpath1+'.so')
        test1 = True
    elif _os.path.exists(libpath1+".dylib"):
        parserlib = _ctypes.cdll.LoadLibrary(libpath1+'.so')
        test1 = True

    #test2 - try /usr/local/lib/libgmadShared.so
    try:
        fulllibname = 'lib'+libname+'.so'
        parserlib   = _ctypes.cdll.LoadLibrary(libpath+fulllibname)
        test2 = True
    except OSError:
        pass

    #test3 - try /usr/local/lib/libgmadShared.dylib
    try:
        fulllibname = 'lib'+libname+'.dylib'
        parserlib   = _ctypes.cdll.LoadLibrary(libpath+fulllibname)
        test3 = True
    except OSError:
        pass

    tests = [test0,test1,test2,test3]
    if tests.count(True) == 0:
        print('LoadLib - cannot find libgmadShared - check paths')
        raise OSError('LoadLib - cannot find libgmadShared - check paths')
    else:
        parserlib.GetNElements.restype   = _ctypes.c_int
        parserlib.GetNElements.argstypes = []
        parserlib.GetType.restype        = _ctypes.c_int
        parserlib.GetType.argtypes       = [_ctypes.c_int]
        parserlib.GetName.restype        = _ctypes.c_char_p
        parserlib.GetName.argtypes       = [_ctypes.c_int]
        parserlib.GetLength.restype      = _ctypes.c_double
        parserlib.GetLength.argtypes     = [_ctypes.c_int]
        parserlib.GetAngle.restype       = _ctypes.c_double
        parserlib.GetAngle.argtypes      = [_ctypes.c_int]
        parserlib.GetKs.restype          = _ctypes.c_double*10
        parserlib.GetKs.argtypes         = [_ctypes.c_int]
        parserlib.GetAper1.restype       = _ctypes.c_double
        parserlib.GetAper1.argtypes      = [_ctypes.c_int]
        parserlib.GetAper2.restype       = _ctypes.c_double
        parserlib.GetAper2.argtypes      = [_ctypes.c_int]
        parserlib.GetAper3.restype       = _ctypes.c_double
        parserlib.GetAper3.argtypes      = [_ctypes.c_int]
        parserlib.GetAper4.restype       = _ctypes.c_double
        parserlib.GetAper4.argtypes      = [_ctypes.c_int]
        parserlib.GetApertureType.restype      = _ctypes.c_char_p
        parserlib.GetApertureType.argtypes     = [_ctypes.c_int]
        parserlib.GetBeampipeThickness.restype  = _ctypes.c_double
        parserlib.GetBeampipeThickness.argtypes = [_ctypes.c_int]

        return parserlib

class Lattice(object):
    """
    BDSIM Gmad parser lattice.

    Use this class to load a bdsim input file using the BDSIM parser (GMAD)
    and then interrogate it.  You can use this to regenerate a lattice with
    less information for example

    >>> a = Lattice("filename.gmad")

    or

    >>> a = Lattice()
    >>> a.Load("filename.gmad")
    >>> a  # this will tell you some basic details
    >>> print(a) # this will print out the full lattice

    """
    def __init__(self,filename=None):
        self._LoadLib()
        if filename != None:
            self.Load(filename)

    def _LoadLib(self):
        self._parserlib = _LoadLib()

    def __repr__(self):
        s =  'BDSIM Gmad Lattice\n'
        s += str(self.nelements) + ' elements in the lattice\n'
        s += str(self.nelements_unique) + ' unique elements\n'
        return s

    def __str__(self):
        s = self._GenerateReprString()
        return s

    def __iter__(self):
        self._iterindex = -1
        if hasattr(self,"lattice") == False:
            raise IOError("No gmad file loaded - cannot iterate")
        return self

    def __len__(self):
        return self.nelements

    def next(self):
        if self._iterindex == len(self.sequence)-1:
            raise StopIteration
        self._iterindex += 1
        return self.lattice[self._iterindex]

    def __getitem__(self,index):
        if type(index) == str:
            ind = self.GetIndexOfElementNamed(index)
            return self.lattice[ind]
        else:
            return self.lattice[index]

    def Load(self, filename):
        """
        Load the BDSIM input file and parse it using the
        BDSIM parser (GMAD).
        """
        if ".gmad" not in filename:
            raise IOError("Not .gmad file - incorrect file type")
        self._parserlib.GmadParser_c(filename)
        self.nelements        = self._parserlib.GetNElements()
        self.sequence         = self.GetAllNames()
        self.names            = list(set(self.sequence))
        self.nelements_unique = len(self.names)
        self.ParseLattice()

    def ParseLattice(self):
        """
        Put lattice data into python data structure
        """
        self.lattice = []
        self.elements = {}
        sposstart = 0.0
        sposend   = 0.0
        for i in range(self.nelements):
            d = self.GetElement(i)
            d['Index']     = i
            d['SPosStart'] =  sposstart
            sposend        += d['Length']
            d['SPosEnd']   =  sposend
            d['SPosMid']   = sposstart + d['Length']/2
            sposstart      += d['Length']
            self.lattice.append(d)

        self.smax = sposend

    def GetName(self,index):
        return str(self._parserlib.GetName(index))

    def GetType(self,index):
        return int(self._parserlib.GetType(index))

    def GetLength(self,index):
        return float(self._parserlib.GetLength(index))

    def GetAngle(self,index):
        return float(self._parserlib.GetAngle(index))

    def GetAper1(self,index):
        return float(self._parserlib.GetAper1(index))

    def GetAper2(self,index):
        return float(self._parserlib.GetAper2(index))

    def GetAper3(self,index):
        return float(self._parserlib.GetAper3(index))

    def GetAper4(self,index):
        return float(self._parserlib.GetAper4(index))

    def GetApertureType(self, index):
        return str(self._parserlib.GetApertureType(index))

    def GetAllNames(self):
        allnames = []
        for i in range(self.nelements):
            allnames.append(self.GetName(i))
        return allnames

    def GetIndexOfElementNamed(self,elementname):
        try:
            ind = self.sequence.index(elementname)
            return ind
        except ValueError("Unknown element name: "+str(elementname)):
            pass

    def GetKs(self,index):
        ks_c = self._parserlib.GetKs(index)
        KsArrayType = _ctypes.c_double*5
        array_pointer = _ctypes.cast(ks_c, _ctypes.POINTER(KsArrayType))
        ks = _np.frombuffer(array_pointer.contents)
        ks[ks < 1e-20] = 0 #replace dbl_min values
        keys = ['ks','k1','k2','k3','k4']
        d = dict(zip(keys,ks))
        return d

    def GetElement(self, i) :
        d = {}
        d['Name']   = self.GetName(i)
        d['Type']   = _GetTypeName(self.GetType(i))
        d['Length'] = self.GetLength(i)
        d['Angle']  = self.GetAngle(i)
        d['Aper1']  = self.GetAper1(i)
        d['Aper2']  = self.GetAper2(i)
        d['Aper3']  = self.GetAper3(i)
        d['Aper4']  = self.GetAper4(i)
        d['ApertureType'] = self.GetApertureType(i)
        d['Ks']     = self.GetKs(i)
        return d

    def _ReprHeader(self):
        h =  '\n'
        h += 'Index'.ljust(7)
        h += 'Name'.ljust(20)
        h += 'Type'.ljust(13)
        h += 'L (m)'.ljust(10)
        h += 'Angle (mrad)'.ljust(15)
        h += 'S Start(m)'.ljust(15)
        h += 'S End  (m)'.ljust(15)
        h += '\n\n'
        return h

    def _GenerateReprString(self,includeheaderlines=True):
        s = ''
        if hasattr(self,'lattice') == False:
            self.ParseLattice()
        for e in self.lattice:
            if (e['Index']%20 == 0) and includeheaderlines:
                s += self._ReprHeader()
            s += str(e['Index']).ljust(7)
            s += str(e['Name']).ljust(20)
            s += str(e['Type']).ljust(13)
            s += str(round(e['Length'],5)).ljust(10)
            s += str(round(e['Angle']*1000.0,9)).ljust(15)
            s += str(round(e['SPosStart'],5)).ljust(15)
            s += str(round(e['SPosEnd']  ,5)).ljust(15)
            s += '\n'
        return s

    def Print(self,includeheaderlines=True):
        s = self._GenerateReprString(includeheaderlines)
        print(s)

    def PrintZeroLength(self,includeheaderlines=True) :
        """
        Print elements with zero length with s location
        """
        if hasattr(self,'lattice') == False:
            self.ParseLattice()
        s = ''
        self.nzerolength = 0
        for e in self.lattice:
            if e['Length'] < 1e-12:
                #it's zero length!
                self.nzerolength += 1
                if (self.nzerolength%20 == 0) and includeheaderlines:
                    s += self._ReprHeader()
                s += str(e['Index']).ljust(7)
                s += str(e['Name']).ljust(20)
                s += str(e['Type']).ljust(13)
                s += str(round(e['Length'],5)).ljust(10)
                s += str(round(e['Angle']*1000.0,9)).ljust(15)
                s += str(round(e['SPosStart'],5)).ljust(15)
                s += str(round(e['SPosEnd']  ,5)).ljust(15)
                s += '\n'
        print(s)

    def IndexFromNearestS(self, S):
        """
        IndexFromNearestS(S)

        return the index of the beamline element clostest to S
        """
        if S < 0 or S > self.smax:
            raise ValueError("S is out of bounds!")
        for index, _ in enumerate(self):
            if S > self[index]["SPosStart"] and S < self[index]["SPosEnd"]:
                return index

    def GetColumn(self, column):
        data = []
        for element in self.lattice:
            data.append(element[column])
        return data

class GmadFile(object) :
    """
    Class to determine parameters and gmad include structure
    """
    def __init__(self, fileName) :

        pass


class GmadFileBeam(object) :
    """
    Class to load a gmad options file to a buffer and modify the contents
    """

    def __init__(self, fileName) :
        pass

class GmadFileOptions(object) :
    """
    Class to load a gmad options file to a buffer and modify the contents
    """

    def __init__(self, fileName) :
        pass


class GmadFileComponents(object) :
    """
    Class to load a gmad components file to a buffer and modify the contents

    Example :
    python> g = pybdsim.Gmad.GmadFileComponents("./atf2_components.gmad")
    python> g.change("KEX1A","l","10")
    python> g.write("./atf2_components.gmad")

    """
    def __init__(self, fileName) :
        '''Load gmad file'''

        # open file and read all of it
        self.fileName = fileName
        f = open(fileName)
        self.s = f.read(-1)
        f.close()

        # regular expressions used for matching

        # For finding elements
        self.elementNameRe = "([a-zA-Z0-9_-]+)\s*:\s*([,a-zA-Z0-9=.\*\s-]+);"
        # For finding a known element
        self.elementRe     = "\s*:\s*([,a-zA-Z0-9=.\*\s-]+);"
        # For extracting the parameters for an element
        self.elementValRe  = "([a-zA-Z0-9_-]+)=([a-zA-Z-0-9.eE\*]+)"
        # For extracting element type
        self.elementNameRe1    = "([a-zA-Z0-9_-]+)\s*:\s*([,a-zA-Z0-9._-]+),([,a-zA-Z0-9=.\*\s-]+);"

        # determine element names in file
        self.elementNames();

    def elementNames(self) :
        '''Make a list of element names, stored in self.elementNameList'''
        self.elementNameList = []
        self.elementTypeList = []
        for l in self.s.split("\n") :
            rem = _re.search(self.elementNameRe1,l)
            try :
                self.elementNameList.append(rem.group(1))
                self.elementTypeList.append(rem.group(2))
            except AttributeError :
                pass

    def findElement(self,elementName) :
        '''Returns the start and end (inclusive location of the element lines as a tuble (start,end)'''

        # See if element is in list
        self.elementNameList.index(elementName)

        # Extract value
        rem = _re.search(elementName+self.elementRe,self.s)

        return [self.s.find(rem.group(0)), len(rem.group(0)), rem.group(0)]

    def parseElement(self,elementString) :
        '''Create element dictionary from element'''
        rem = _re.findall(self.elementValRe,elementString)

        d = {}
        for v in rem :
            d[v[0]] = v[1]

        return d

    def change(self,element,parameter, value) :
        '''Edit element dictionary'''
        e = self.findElement(element)
        d = self.parseElement(e[2])

        d[parameter] = value

        # form new string
        s = e[2][0:e[2].find(",")]
        for k in d.keys() :
            s+=", "+k+"="+d[k]
        s+=";"

        # replace string of gmad file
        self.s = self.s.replace(e[2],s)

    def getParameter(self,element,parameter) :
        '''Edit element dictionary'''
        e = self.findElement(element)
        d = self.parseElement(e[2])
        return d[parameter]

    def getType(self,element) :
        i = self.elementNameList.index(element)
        return self.elementTypeList[i]

    def write(self,fileName) :
        f = open(fileName,"w+")
        f.write(self.s)
        f.close()

def _GetTypeName(typeenum):
    """
    Convert the enum to a proper name"
    """
    try:
        return _BDSIMTypes[typeenum]
    except KeyError:
        return 'Unknown'

_BDSIMTypes = {
    -1  : 'None',
    1   : 'Marker',
    2   : 'Drift',
    63  : 'PCL Drift',
    3   : 'Drift',
    4   : 'SBend',
    5   : 'Quadrupole',
    6   : 'Sextupole',
    7   : 'Octupole',
    9   : 'Multipole',
    10  : 'Solenoid',
    11  : 'Line',
    -11 : 'Reversed Line',
    12  : 'Collimator',
    13  : 'E-Collimator',
    62  : 'Mu-Spoiler',
    14  : 'R-Collimator',
    15  : 'Laser',
    16  : 'Material',
    17  : 'RBend',
    18  : 'Atom',
    19  : 'Sequence',
    21  : 'Screen',
    22  : 'Awake Screen',
    31  : 'VKick',
    32  : 'HKick',
    41  : 'Sampler',
    42  : 'CSampler',
    43  : 'Dump',
    51  : 'Gas',
    52  : 'Tunnel',
    61  : 'Transform3D',
    98  : 'Teleporter',
    99  : 'Terminator'
    }
