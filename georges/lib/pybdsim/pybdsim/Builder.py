# pybdsim.Builder - tools to build bdsim lattices
# Version 1.0
# L. Nevay
# laurie.nevay@rhul.ac.uk

"""
Builder

Build generic machines for bdsim. You can create a lattice
using one of the predefined simple lattices or by adding
many pieces together of your own design. Finally, output
the gmad files required.

Classes:
Element - beam line element that always has name,type and length
Machine - a list of elements

"""
from . import Beam as _Beam
from . import Options as _Options
from . import Writer as _Writer
from . import _General
from ._General import IsFloat as _IsFloat
import math as _math
import time as _time
import os as _os
import numpy as _np
import copy as _copy
import textwrap as _textwrap

bdsimcategories = [
    'marker',
    'drift',
    'rbend',
    'sbend',
    'quadrupole',
    'sextupole',
    'octupole',
    'decapole',
    'multipole',
    'thinmultipole',
    'rfcavity',
    'rcol',
    'ecol',
    'muspoiler',
    'solenoid',
    'hkicker',
    'vkicker',
    'kicker',
    'tkicker',
    'transform3d',
    'element',
    'line',
    'matdef',
    'laser',
    'gas',
    'spec',
    'degrader',
    'shield',
    'gap',
    'placement',
    ]

class ElementBase(dict):
    """
    A class that represents an element / item in an accelerator beamline.
    Printing or string conversion produces the BDSIM syntax.

    This class provides the basic dict(ionary) inheritance and functionality 
    and the representation that allows modification of existing parameters
    of an already declared item.

    """
    def __init__(self, name, isMultipole=False, **kwargs):
        dict.__init__(self)
        self['name']      = name
        self.name         = name
        self._isMultipole = isMultipole
        self._keysextra   = []
        self.__Update(kwargs)

    def __Update(self,d):
        dict.update(d)
        for key,value in d.items():
            if value == "" :
                continue
            if type(value) == tuple and self._isMultipole:
                self[key] = value
            elif type(value) == tuple:
                #use a tuple for (value,units)
                self[key] = (float(value[0]),value[1])
            elif isinstance(value, (float, int)):
                #must be a number
                if 'aper' in str.lower(key) and value < 1e-6:
                    continue
                else:
                    self[key] = value
            else:
                #must be a string
                self[key] = '"'+value+'"'
            self._keysextra.append(str(key)) #order preserving

    def keysextra(self):
        #so behaviour is similar to dict.keys()
        return self._keysextra

    def __repr__(self):
        s = ''
        s += self.name + ': '
        for i,key in enumerate(self._keysextra):
            if i > 0:
                s += ", "
            if type(self[key]) == tuple and self._isMultipole:
                s += key + '=' + '{'+(','.join([str(s) for s in self[key]]))+'}'
            elif type(self[key]) == tuple:
                s += key + '=' + str(self[key][0]) + '*' + str(self[key][1])
            else:
                s += key + '=' + str(self[key])
        s += ';\n'
        return s

class Element(ElementBase):
    """
    Element - an element / item in an accelerator beamline. Very similar to a
    python dict(ionary) and has the advantage that built in printing or string
    conversion provides BDSIM syntax.

    Element(name,type,**kwargs)

    >>> a = Element("d1", "drift", l=1.3)
    >>> b = Element("qx1f", "quadrupole", l=(0.4,'m'), k1=0.2, aper1=(0.223,'m'))
    >>> print(b)
    qx1f: quadrupole, k1=0.2, l=0.4*m, aper1=0.223*m;
    >>> str(c)
    qx1f: quadrupole, k1=0.2, l=0.4*m, aper1=0.223*m\\n;
    
    A beam line element must ALWAYs have a name, and type.
    The keyword arguments are specific to the type and are up to
    the user to specify - these should match BDSIM GMAD syntax.

    The value can be either a single string or number or a python tuple
    where the second entry must be a string (shown in second example).
    Without specified units, the parser assumes S.I. units.

    An element may also be multiplied or divided.  This will scale the 
    length and angle appropriately.

    >>> c = Element('sb1', 'sbend', l=(0.4,'m'), angle=0.2)
    >>> d = c/2
    >>> print(d)
    sb1: sbend, l=0.2*m, angle=0.1;

    This inherits and extends ElementBase that provides the basic dictionary
    capabilities.  It adds the requirement of type / category (because 'type' is
    a protected keyword in python) as well as checking for valid BDSIM types.
    """
    def __init__(self, name, category, **kwargs):
        if category not in bdsimcategories:
            raise ValueError("Not a valid BDSIM element type")
        if (category == 'thinmultipole') or (category == 'multipole'):
            ElementBase.__init__(self,name,isMultipole=True,**kwargs)
        else:
            ElementBase.__init__(self,name,**kwargs)
        self['category'] = category
        self.category    = category
        self.length      = 0.0 #for book keeping only
        self.__Update()
        self._UpdateLength()
        
    def __Update(self, d=None):
        if d != None:
            ElementModifier.__update(self,d)

    def _UpdateLength(self):
        if 'l' in self:
            if type(self['l']) == tuple:
                ll = self['l'][0]
            else:
                ll = self['l']
            self.length += float(ll)

    def __repr__(self):
        s = "{}: {}".format(self.name, self.category)
        if self._keysextra:
            for key in sorted(self._keysextra):
                if (type(self[key]) == tuple
                    and (self.category != 'thinmultipole')
                    and (self.category != 'multipole')):
                    s += (', ' + key + '=' + str(self[key][0])
                          + '*' + str(self[key][1]))
                elif (type(self[key]) == tuple
                      and ((self.category == 'thinmultipole')
                           or (self.category == 'multipole'))):
                    s += (', ' + key + '=' + '{'
                          + (','.join([str(s) for s in self[key]]))+'}')
                else:
                    s += ", {}={}".format(key, self[key])
        s = '\n\t'.join(_textwrap.wrap(s))
        s += ';\n'
        return s

    def __mul__(self, factor):
        newElement = _copy.copy(self)
        newElement.length *= factor
        newElement['l'] = factor * float(self['l'])
        if 'angle' in self:
            newElement['angle'] = factor * float(self['angle'])
        return newElement

    def __div__(self, factor):
        return self.__mul__(float(1./factor))

class ElementModifier(ElementBase):
    """
    A class  to MODIFY an already defined element in a gmad file by appending an
    updated definition. Using this alone in BDSIM will result in an 
    undefined type error. This class is particularly useful for creating
    a strength file.

    # define an element
    >>> a = Element('qf1', 'quadrupole', l=0.3, k1=0.00345)
    >>> b = ElementModifier('qf1',k1=0.0245)
    >>> f = open('mylattice.gmad', 'w')
    >>> f.write(str(a))
    >>> f.write(str(b))
    >>> f.close()

    cat mylattice.gmad
    qf1, quadrupole, l=0.3, k1=0.00345;
    qf1, k1=0.0245

    This results in the quadrupole strength k1 in this example being 
    changed to 0.0245.
    """
    def __init__(self, name, isMultipole=False, **kwargs):
        if len(kwargs) == 0:
            raise ValueError("Error: must specify at least 1 keyword argument")
        ElementBase.__init__(self, name, isMultipole, **kwargs)

class Line(list):
    """
    A class that represents a :class:`list` of :class:`Elements`

    Provides ability to print out the sequence or define all 
    the components.

    Example:

    >>> d1 = Element("drift1", "drift", l=1.3)
    >>> q1 = Element("q1", "quadrupole", l=0.4, k1=4.5)
    >>> a = Line([d1,q1])

    """
    def __init__(self,name,*args):
        for item in args[0]:
            if type(item) != Element:
                raise TypeError("Line is a list of Elements")
        list.__init__(self,*args)
        self.name   = name
        self.length = 0.0 
        for item in args[0]:
            self.length += item.length
        
    def __repr__(self):
        s = ''
        for item in self:
            s += str(item)+'\n' #uses elements __repr__ function
        s += self.name+ ': line=('

        s += ', '.join([item.name for item in self]) + ');'
        s = '\n\t'.join(_textwrap.wrap(s))
        s += "\n"
        return s

    def DefineConstituentElements(self):
        """
        Return a string that contains the lines required
        to define each element in the :class:`Line`.

        Example using predefined Elements name 'd1' and 'q1':

        >>> l = Line([d1,q1])
        >>> f = open("file.txt", "w")
        >>> f.write(DefineConsituentElements())
        >>> f.write(l)
        >>> f.close()
        """
        s = ''
        for item in self:
            s += str(item) #uses elements __repr__ function
        return s

class ApertureModel(dict):
    """
    A class that produces the aperture representation of an element. Only non-zero
    values are written for the aperture parameters. Includes parameter checking.
    """
    def __init__(self, apertureType='circular', aper1=0.1, aper2=0, aper3=0, aper4=0):
        dict.__init__(self)
        allowedTypes = [
            'circular',
            'elliptical',
            'rectangular',
            'lhc',
            'lhcdetailed',
            'rectellipse'
        ] # maintain order for tests further down!
        madxTypes = {
            'circle'      : 'circular',
            'ellipse'     : 'elliptical',
            'rectangle'   : 'rectangular',
            'lhcscreen'   : 'lhcdetailed',
            'marguerite'  : None,
            'rectellipse' : 'rectellipse',
            'racetrack'   : None,
            'filename'    : None
            }
        atL = str.lower(apertureType)
        if atL not in allowedTypes and atL not in madxTypes.keys():
            print('Allowed aperture types are: ', allowedTypes, madxTypes.keys())
            raise ValueError("Invalid aperture type: "+str(apertureType))

        if atL in madxTypes.keys():
            self['apertureType'] = madxTypes[atL]
            if self['apertureType'] == None:
                print('Unsupported type: :',self.apertureType,'" - replacing with elliptical')
                self['apertureType'] = 'elliptical'
        else:
            self['apertureType'] = apertureType

        if self['apertureType'] in allowedTypes[1:] and aper2 == 0:
            print('For aperture type "',self['apertureType'],'" at least aper1 and aper2 must be specified')
            raise ValueError("Too few aperture parameters supplied")
        
        self['aper1'] = aper1
        self['aper2'] = aper2
        self['aper3'] = aper3
        self['aper4'] = aper4

    def __repr__(self):
        # aper1 is always present at least.
        out = ('apertureType="{}*m", aper1={}*m').format(self["apertureType"],
                                                         self["aper1"])
        # Append any non-zero apertures.
        for i in [2,3,4]:
            aperKey = "aper{}".format(i)
            aperValue = self[aperKey]
            if aperValue > 0:
                out += ", {}={}*m".format(aperKey, aperValue)
        return out

class Sampler:
    """
    A sampler is unique in that it does not have a length unlike every
    :class:`Element` hence it needs its own class to produce its 
    representation.
    """
    def __init__(self,name):
        self.name = name

    def __repr__(self):
        if self.name == 'all':
            return 'sample, all;\n'
        else:
            return 'sample, range='+self.name+';\n'

class Machine:
    """
    A class represents an accelerator lattice as a sequence of
    components. Member functions allow various lattice components
    to be append to the sequence of the machine. This class allows
    the user to programatically create a lattice and write the 
    BDSIM gmad representation of it.

    Example:
    
    >>> a = Machine()
    >>> a.AddDrift('mydrift', l=1.3)
    >>> a.Write("lattice.gmad")

    Example with Sychrotron rescaling:
    
    >>> a = Machine(sr=True, energy0=250,charge=-1)
    >>> a.AddDipole('sb1','sbend',length=1.0,1e-5)
    >>> a.AddDrift('dr1',length=1)
    >>> a.AddDipole('sb2','sbend',length=1.0,1e-5)
    >>> a.AddDrift("dr2",length=1)

    Caution: adding an element of the same name twice will result the 
    element being added only to the sequence again and not being
    redefined - irrespective of if the parameters are different. If
    verbose is used (True), then a warning will be issued.    
    
    """
    def __init__(self,verbose=False, sr=False, energy0=0.0, charge=-1.0):
        self.verbose   = verbose
        self.sequence  = []
        self.elements  = []
        self.elementsd = {}
        self.samplers  = []
        self.length    = 0.0
        self.angint    = 0.0
        self.bias      = []
        self.beam      = _Beam.Beam()
        self.options   = None
        self.energy0   = energy0
        self.energy    = []
        self.lenint    = []
        self.sr        = sr
        self.energy.append(energy0)
        self.charge = charge

    def __repr__(self):
        s = ''
        s += 'pybdsim.Builder.Machine instance\n'
        s += str(len(self.sequence)) + ' items in sequence\n'
        s += str(len(self.elements)) + ' unique elements defined\n'
        return s

    def __iter__(self):
        self._iterindex = -1
        return self

    def next(self):
        if self._iterindex == len(self.sequence)-1:
            raise StopIteration
        self._iterindex += 1
        return self.elementsd[self.sequence[self._iterindex]]
        
    def __getitem__(self,name):
        if _IsFloat(name):
            return self.elementsd[self.sequence[name]]
        else:
            return self.elementsd[name]

    def __len__(self):
        print('Number of unique elements:      ',len(self.elementsd.keys()))
        print('Number of elements in sequence: ',len(self.sequence),' <- returning this')
        return len(self.sequence)

    def GetIntegratedAngle(self):
        """
        Get the cumulative angle of all the bends in the machine. This is therefore the difference
        in angle between the entrance and exit vectors. All angles are assumed to be in the horizontal
        plane so this will not be correct for rotated dipoles.
        """
        return self.angint

    def GetIntegratedLength(self):
        """
        Get the integrated length of all the components.
        """
        return self.length

    def Append(self, object):
        if type(object) not in (Element,Line):
            raise TypeError("Only Elements or Lines can be added to the machine")
        elif object.name not in self.elementsd.keys():
            #hasn't been used before - define it
            if type(object) is Line:
                for element in object:
                    self.Append(object)
            else:
                self.elements.append(object)
                self.elementsd[object.name] = object
        else:
            if self.verbose:
                print("Element of name: ",object.name," already defined, simply adding to sequence")
        #finally add it to the sequence
        if object.category is not "placement":
            self.sequence.append(object.name)

        self.length += object.length
        self.lenint.append(self.length)

        # list of elements that produce SR
        elementsSR = ["sbend", "rbend"]

        # update energy if correct element category and has finite length.
        if object.category in elementsSR and object.length > 0:
            if 'angle' in object:
                ang = object['angle']
                if type(ang) == tuple:
                    ang = ang[0]
                else:
                    ang = ang
            elif 'B' in object:
                # Assume a beam instance has been added to machine...
                if (self.beam['particle'] == "e-") or (self.beam['particle'] == "e+"):
                    pMass = 0.000511
                elif (self.beam['particle'] == "proton"):
                    pMass = 0.938
                else:
                    pMass = 0
                if (self.energy[-1] > 0):
                    brho = 3.3356 * _np.sqrt(self.energy[-1]**2 - pMass**2)
                    ang  = object['B'] * object.length / brho
                else:
                    ang  = 0
            else:
                ang = 0

            self.angint += ang

            # bend so SR generated. Recompute beam energy
            energy = self.energy[-1]
            self.energy.append(energy-14.1e-6*ang**2/object.length*energy**4)
        else :
            self.energy.append(self.energy[-1])

    def SynchrotronRadiationRescale(self):
        """
        Rescale all component strengths for SR
        """
        ielement = 1
        for element in self.elements:
            # energyave = (self.energy[ielement]+self.energy[ielement-1])/2.0
            energyave = self.energy[ielement]
            # print energyave
            if element.category == 'rbend' or element.category == 'sbend' :
                angle  = element['angle']
                length = element['l']

                # insert magnetic field value after angle
                element._keysextra.insert(element._keysextra.index('angle')+1,'B')
                # consistent calculation with BDSIM
                element['B'] = self.charge*energyave/0.299792458*angle/length
            elif element.category == 'quadrupole' :
                element['k1'] = energyave / self.energy0 * element['k1']
            elif element.category == 'sextupole' :
                element['k2'] = energyave / self.energy0 * element['k2']
            elif element.category == 'octupole':
                element['k3'] = energyave / self.energy0 * element['k3']
            elif element.category == 'decupole':
                element['k4'] = energyave / self.energy0 * element['k4']
            elif element.category == 'multipole' :
                pass
            ielement += 1

    def Write(self, filename, verbose=False, overwrite=True):
        """
        Write the machine to a series of gmad files.

        kwargs:
        overwrite : Do not append an integer to the basefilename if
        already exists, instead overwrite existing files.

        """
        if self.sr :
            self.SynchrotronRadiationRescale()

        verboseresult = verbose or self.verbose
        writer = _Writer.Writer()
        writer.WriteMachine(self,filename,verboseresult)

    def AddBias(self, biasobject):
        self.bias.append(biasobject)

    def AddBeam(self, beam=None):
        """
        Assign a beam instance to this machine. If no Beam instance is provided,
        a reference distribution is used.
        """
        if beam == None:
            self.beam = _Beam.Beam()
        elif type(beam) != _Beam.Beam:
            raise TypeError("Incorrect type - please provide pybdsim.Beam.Beam instance")
        self.beam = beam

    def AddOptions(self, options=None):
        """
        Assign an options instance to this machine. 
        """
        if type(options) != _Options.Options:
            raise TypeError("Incorrect type - please provide pybdsim.Options.Options instance")
        self.options = options
        
    def AddMarker(self, name='mk'):
        """ 
        Add a marker to the beam line.
        """
        if self.verbose:
            print('AddMarker> ',name)
        self.Append(Element(name,'marker'))
    
    def AddDrift(self, name='dr', length=0.1, **kwargs):
        """
        Add a drift to the beam line
        """
        if self.verbose:
            print('AddDrift>  ',name,' ',length,' ',kwargs)
        if length < 1e-12:
            self.AddMarker(name)
        else:
            self.Append(Element(name,'drift',l=length,**kwargs))

    def AddDipole(self, name='dp', category='sbend', length=0.1,
                  angle=None, b=None, **kwargs):
        """
        AddDipole(category='sbend')

        category - 'sbend' or 'rbend' - sector or rectangular bend

        """
        if category not in ['sbend','rbend']:
            s = 'Invalid category ' + category
            raise ValueError(s)
        if angle is None and b is None:
            raise TypeError('angle or b must be specified for a dipole')
        elif angle != None:
            self.Append(Element(name,category,l=length,angle=angle,**kwargs))
        else:
            self.Append(Element(name,category,l=length,B=b,**kwargs))

    def AddQuadrupole(self, name='qd', length=0.1, k1=0.0, **kwargs):
        self.Append(Element(name,'quadrupole',l=length,k1=k1,**kwargs))
        
    def AddSextupole(self, name='sx', length=0.1, k2=0.0, **kwargs):
        self.Append(Element(name,'sextupole',l=length,k2=k2,**kwargs))

    def AddOctupole(self, name='oc', length=0.1, k3=0.0, **kwargs):
        self.Append(Element(name,'octupole',l=length,k3=k3,**kwargs))

    def AddDecapole(self, name='dc', length=0.1, k4=0.0, **kwargs):
        self.Append(Element(name,'decapole',l=length,k4=k4,**kwargs))

    def AddMultipole(self, name='mp', length=0.1, knl=(0), ksl=(0), tilt=0.0, **kwargs):
        if length > 1e-12:
            self.Append(Element(name,'multipole',l=length, knl=knl, ksl=ksl, tilt=tilt, **kwargs))
        else :
            self.AddMarker(name)

    def AddThinMultipole(self, name='mp', knl=(0), ksl=(0), tilt=0.0, **kwargs):
        self.Append(Element(name,'thinmultipole', knl=knl, ksl=ksl, tilt=tilt, **kwargs))

    def AddRFCavity(self, name='arreff', length=0.1, gradient=10, **kwargs) :
        self.Append(Element(name,'rfcavity',l=length, gradient=gradient, **kwargs))
        
    def AddRCol(self, name='rc', length=0.1, xsize=0.1, ysize=0.1, **kwargs):
        d = {}
        for k,v in kwargs.items():
            if 'aper' not in str(k).lower():
                d[k] = v
        self.Append(Element(name,'rcol',l=length,xsize=xsize,ysize=ysize,**d))

    def AddDegrader(self, length=0.1, name='deg', nWedges=1, wedgeLength=0.1, degHeight=0.1, materialThickness=None, degraderOffset=None, **kwargs):
        if (materialThickness==None) and (degraderOffset==None):
            raise TypeError('materialThickness or degraderOffset must be specified for a degrader')
        elif materialThickness != None:
            self.Append(Element(name,'degrader',l=length,numberWedges=nWedges,wedgeLength=wedgeLength,
                                degraderHeight=degHeight, materialThickness=materialThickness,**kwargs))
        else:
            self.Append(Element(name,'degrader',l=length,numberWedges=nWedges,wedgeLength=wedgeLength,
                                degraderHeight=degHeight, degraderOffset=degraderOffset,**kwargs))

    def AddMuSpoiler(self, name='mu', length=0.1, b=0.0, **kwargs):
        self.Append(Element(name,'muspoiler',l=length,B=b,**kwargs))

    def AddSolenoid(self, name='sl', length=0.1, ks=0.0, **kwargs):
        self.Append(Element(name,'solenoid',l=length,ks=ks,**kwargs))

    def AddShield(self, name='sh', length=0.1, **kwargs):
        self.Append(Element(name,'shield',l=length,**kwargs))

    def AddLaser(self, length=0.1, name='lsr', x=1, y=0, z=0, waveLength=532e-9, **kwargs):
        self.Append(Element(name,'laser',l=length,x=x,y=y,z=z,waveLength=waveLength,**kwargs))

    def AddTransform3D(self, name='t3d',**kwargs):
        if len(kwargs.keys()) == 0:
            pass
        else:
            self.Append(Element(name,'transform3d',**kwargs))

    def AddECol(self, name='ec', length=0.1, xsize=0.1, ysize=0.1, **kwargs):
        d = {}
        for k,v in kwargs.items():
            if 'aper' not in str(k).lower():
                d[k] = v
        self.Append(Element(name,'ecol',l=length,xsize=xsize,ysize=ysize,**d))
        
    def AddHKicker(self, name='hk', hkick=0.0, **kwargs):
        self.Append(Element(name,'hkicker', hkick=hkick, **kwargs))

    def AddVKicker(self, name='vk', vkick=0.0, **kwargs):
        self.Append(Element(name,'vkicker', vkick=vkick, **kwargs))

    def AddKicker(self, name='kk', hkick=0.0, vkick=0.0, **kwargs):
        self.Append(Element(name,'kicker', hkick=hkick, vkick=hkick, **kwargs))

    def AddTKicker(self, name='tk', hkick=0.0, vkick=0.0, **kwargs):
        self.Append(Element(name,'tkicker', hkick=hkick, vkick=vkick, **kwargs))

    def AddElement(self, name='el', length=0.1, outerDiameter=1, geometryFile="geometry.gdml", **kwargs):
        self.Append(Element(name, 'element',l=length,outerDiameter=outerDiameter,geometryFile=geometryFile, **kwargs))

    def AddPlacement(self, name='pl', geometryFile="gdml:geometry.gdml",
                     phi=0.0, psi=0.0, theta=0.0,
                     x=0.0, y=0.0, z=0.0,
                     **kwargs):
        self.Append(Element(name, 'placement', geometryFile="gdml:" + geometryFile,
                            psi=psi, phi=phi, theta=theta,
                            x=x, y=y, z=z, **kwargs))

    def AddGap(self, name='gp', length=1.0, **kwargs):
        self.Append(Element(name, 'gap', l=length, **kwargs))

    def AddFodoCell(self, basename='fodo', magnetlength=1.0, driftlength=4.0,kabs=0.2,**kwargs):
        """
        AddFodoCell(basename,magnetlength,driftlength,kabs,**kwargs)
        basename     - the basename for the fodo cell beam line elements
        magnetlength - length of magnets in metres
        driftlength  - length of drift segment in metres
        kabs         - the absolute value of the quadrupole strength - alternates between magnets

        \*\*kwargs are other parameters for bdsim - ie material='Fe'
        """
        names = [basename+extrabit for extrabit in ['_qfa','_dra','_qda','_drb','_qfb']]
        items = (
            Element(names[0],'quadrupole',l=magnetlength/2.0,k1=kabs,**kwargs),
            Element(names[1],'drift',l=driftlength),
            Element(names[2],'quadrupole',l=magnetlength,k1=-1.0*kabs,**kwargs),
            Element(names[3],'drift',l=driftlength),
            Element(names[4],'quadrupole',l=magnetlength/2.0,k1=kabs,**kwargs)
            )
        self.Append(Line(basename,items))

    def AddFodoCellSplitDrift(self, basename='fodo', magnetlength=1.0, driftlength=4.0, kabs=0.2,nsplits=10, **kwargs):
        """
        AddFodoCellSplitDrift(basename,magnetlength,driftlength,kabs,nsplits,**kwargs)
        basename - the basename for the fodo cell beam line elements
        magnetlength - length of magnets in metres
        driftlength  - length of drift segment in metres
        kabs         - the absolute value of the quadrupole strength - alternates between magnets
        nsplits      - number of segments drift length is split into 

        Will add qf quadrupole of strength +kabs, then drift of l=driftlength split into 
        nsplit segments followed by a qd quadrupole of strength -kabs and the same pattern
        of drift segments.
        
        nsplits will be cast to an even integer for symmetry purposes.

        \*\*kwargs are other parameters for bdsim - ie aper=0.2
        """
        nsplits = _General.NearestEvenInteger(nsplits)
        splitdriftlength = driftlength / float(nsplits)
        maxn    = int(len(str(nsplits)))
        self.Append(Element(basename+'_qfa','quadrupole',l=magnetlength/2.0,k1=kabs,**kwargs))
        for i in range(nsplits):
            self.Append(Element(basename+'_d'+str(i).zfill(maxn),'drift',l=splitdriftlength))
        self.Append(Element(basename+'_qd','quadrupole',l=magnetlength,k1=-1.0*kabs,**kwargs))
        for i in range(nsplits):
            self.Append(Element(basename+'_d'+str(i).zfill(maxn),'drift',l=splitdriftlength))
        self.Append(Element(basename+'_qfb','quadrupole',l=magnetlength/2.0,k1=kabs,**kwargs))

    def AddFodoCellMultiple(self, basename='fodo', magnetlength=1.0, driftlength=4.0, kabs=0.2, ncells=2, **kwargs):
        ncells = int(ncells)
        maxn   = int(len(str(ncells)))
        for i in range(ncells):
            cellname = basename+'_'+str(i).zfill(maxn)
            self.AddFodoCell(cellname,magnetlength,driftlength,kabs,**kwargs)

    def AddFodoCellSplitDriftMultiple(self, basename='fodo', magnetlength=1.0, driftlength=4.0, kabs=0.2, nsplits=10, ncells=2, **kwargs):
        ncells = int(ncells)
        maxn   = int(len(str(ncells)))
        for i in range(ncells):
            cellname = basename+'_'+str(i).zfill(maxn)
            self.AddFodoCellSplitDrift(cellname,magnetlength,driftlength,kabs,nsplits=10,**kwargs)
            
    def AddSampler(self,*elementnames):
        if elementnames[0] == 'all':
            self.samplers.append(Sampler('all'))
        elif elementnames[0] == 'first':
            self.samplers.append(Sampler(self.elements[0].name))
        elif elementnames[0] == 'last':
            self.samplers.append(Sampler(self.elements[-1].name))
        else:
            for element in elementnames:
                if element not in self.sequence:
                    raise ValueError(element+" is not a valid element in this machine")
                else:
                    self.samplers.append(Sampler(element))


# General scripts below this point

def PrepareApertureModel(rowDictionary, default='circular'):
    rd = rowDictionary # shortcut
    a1 = rd['APER_1']
    a2 = rd['APER_2']
    a3 = rd['APER_3']
    a4 = rd['APER_4']
    if 'APERTYPE' in rd.keys():
        aType = str.lower(rd['APERTYPE'])
    else:
        # no type given - let's guess :(
        if a1 == a2 == a3 == a4:
            aType = 'circular'
            a2 = a3 = a4 = 0 # set 0 to save writing needlessly
        else:
            aType = default
    a = ApertureModel(aType,a1,a2,a3,a4)
    return a

def CreateDipoleRing(filename, ndipoles=60, circumference=100.0, samplers='first'):
    """
    Create a ring composed of only dipoles
    filename
    ncells        - number of cells, each containing 1 dipole and a drift
    circumference - in metres
    samplers      - 'first', 'last' or 'all'
    
    """
    ndipoles = int(ndipoles)
    a            = Machine()
    dangle       = 2.0*_math.pi / float(ndipoles)
    dipolelength = circumference / float(ndipoles)
    for i in range(ndipoles):
        a.AddDipole(length=dipolelength, angle=dangle)
    a.AddSampler(samplers)
    a.Write(filename)

def CreateDipoleDriftRing(filename, ncells=60, circumference=100.0, driftfraction=0.1, samplers='first'):
    """
    Create a ring composed of dipoles and drifts
    filename
    ncells        - number of cells, each containing 1 dipole and a drift
    circumference - in metres
    driftfraction - the fraction of drift in each cell (0.0 < driftfraction < 1.0)
    samplers      - 'first', 'last' or 'all'
    
    """
    ncells = int(ncells)
    if driftfraction > 1.0:
        raise Warning("Fraction of drift must be less than 1.0 -> setting to 0.9")
        driftfraction = 0.9
    if driftfraction < 0.0:
        raise Warning("Fraction of drift must be greater than 1.0 -> setting to 0.1")
        driftfraction = 0.1
    a            = Machine()
    dangle       = 2.0*_math.pi / float(ncells)
    clength      = circumference / float(ncells)
    driftlength  = clength * driftfraction
    dipolelength = clength - driftlength
    for i in range(ncells):
        a.AddDipole(length=dipolelength*0.5, angle=dangle*0.5)
        a.AddDrift(length=driftlength)
        a.AddDipole(length=dipolelength*0.5, angle=dangle*0.5)
    a.AddSampler(samplers)
    a.Write(filename)

def CreateDipoleFodoRing(filename, ncells=60, circumference=200.0, samplers='first'):
    """
    Create a ring composed of fodo cells with 2 dipoles per fodo cell.

    filename
    ncells         - number of fodo+dipole cells to create
    circumference  - circumference of machine in metres
    samplers       - 'first','last' or 'all'
    
    Hard coded to produce the following cell fractions:
    50% dipoles
    20% quadrupoles
    30% beam pipe / drift
    """
    a       = Machine()
    cangle  = 2.0*_math.pi / ncells
    clength = float(circumference) / ncells
    #dipole = 0.5 of cell, quads=0.2, drift=0.3, two dipoles
    #dipole:
    dl  = clength * 0.25
    da  = cangle * 0.5
    #quadrupole:
    ql  = clength * 0.2 * 0.5
    k1  = SuggestFodoK(ql,dl)
    #drift:
    drl = clength * 0.3 * 0.25
    #naming
    nplaces  = len(str(ncells))
    basename = 'dfodo_'
    for i in range(ncells):
        cellname = basename + str(i).zfill(nplaces)
        a.AddQuadrupole(cellname+'_qd_a',ql*0.5,k1)
        a.AddDrift(cellname+'_dr_a',drl)
        a.AddDipole(cellname+'_dp_a','sbend',dl,da)
        a.AddDrift(cellname+'_dr_b',drl)
        a.AddQuadrupole(cellname+'_qf_b',ql,k1*-1.0)
        a.AddDrift(cellname+'_dr_c',drl)
        a.AddDipole(cellname+'_dp_b','sbend',dl,da)
        a.AddDrift(cellname+'_dr_d',drl)
        a.AddQuadrupole(cellname+'_qd_c',ql*0.5,k1)
    a.AddSampler(samplers)
    a.Write(filename)
    
def CreateFodoLine(filename, ncells=10, driftlength=4.0, magnetlength=1.0, samplers='all',**kwargs):
    """
    Create a FODO lattice with ncells.
    
    ncells       - number of fodo cells
    driftlength  - length of drift segment in between magnets
    magnetlength - length of quadrupoles
    samplers     - 'all','first' or 'last'
    \*\*kwargs   - kwargs to supply to quadrupole constructor
    
    """
    ncells = int(ncells)
    a      = Machine()
    k1     = SuggestFodoK(magnetlength,driftlength)
    a.AddFodoCellSplitDriftMultiple(magnetlength=magnetlength,driftlength=driftlength,kabs=k1,nsplits=10,ncells=ncells,**kwargs)
    a.AddSampler(samplers)
    a.Write(filename)

def SuggestFodoK(magnetlength,driftlength):
    """
    SuggestFodoK(magnetlength,driftlength)

    returns k1 (float) value for matching into next quad in a FODO cell.
    f = 1/(k1 * magnetlength) = driftlength -> solve for k1

    Note the convention in pybdsim.Builder is that the quadrupoles in
    the fodo cell are split in two.  So this is in fact half the integrated
    k you need.  This matches with the other functions in Builder.

    """
    return 1.0 / (float(magnetlength)*(float(magnetlength) + float(driftlength)))

def WriteMachine(machine, filename, verbose=False):
    """
    WriteMachine(machine(machine),filename(string),verbose(bool))

    Write a machine to disk. This writes several files to make the
    machine, namely:
    
    filename_components.gmad - component files (max 10k per file)
    filename_sequence.gmad   - lattice definition
    filename_samplers.gmad   - sampler definitions (max 10k per file)
    filename_options.gmad    - options (TO BE IMPLEMENTED)
    filename.gmad            - suitable main file with all sub 
                               files in correct order
    
    these are prefixed with the specified filename / path
    
    """
    
    if not isinstance(machine,Machine):
        raise TypeError("Not machine instance")
    
    elementsperline = 100 #number of machine elements per bdsim line (not text line)
    
    #check filename
    if filename[-5:] != '.gmad':
        filename += '.gmad'

    #check for directory and make it if not:
    if '/' in filename:
        directory = '/'.join(filename.split('/')[:-1]) #strip the filename off
        if not _os.path.exists(directory):
            _os.system("mkdir -p " + directory)
    
    #check if file already exists
    ofilename = filename
    filename = _General.GenUniqueFilename(filename)
    if filename != ofilename:
        print('Warning, chosen filename already exists - using filename: ',filename.split('.')[0])
    basefilename = filename[:-5] #everything before '.gmad'

    #prepare names
    files         = []
    fn_main       = basefilename + '.gmad'
    fn_components = basefilename + '_components.gmad'
    fn_sequence   = basefilename + '_sequence.gmad'
    fn_samplers   = basefilename + '_samplers.gmad'
    fn_beam       = basefilename + '_beam.gmad'
    fn_options    = basefilename + '_options.gmad'
    fn_bias       = basefilename + '_bias.gmad'
    timestring = '! ' + _time.strftime("%a, %d %b %Y %H:%M:%S +0000", _time.gmtime()) + '\n'

    #write bias if it exists
    if len(machine.bias) > 0:
        f = open(fn_bias,'w')
        files.append(fn_bias)
        for bias in machine.bias:
            f.write(str(bias))
        f.close()
    
    #write component files
    f = open(fn_components, 'w')
    files.append(fn_components)
    f.write(timestring)
    f.write('! pybdsim.Builder Lattice \n')
    f.write('! COMPONENT DEFINITION\n\n')
    for element in machine.elements:
        f.write(str(element))
    f.close()

    #write lattice sequence
    f = open(fn_sequence,'w')
    files.append(fn_sequence)
    f.write(timestring)
    f.write('! pybdsim.Builder Lattice \n')
    f.write('! LATTICE SEQUENCE DEFINITION\n\n')
    sequencechunks = _General.Chunks(machine.sequence,elementsperline)
    linelist = []
    ti = 0
    for line in _General.Chunks(machine.sequence,elementsperline):
        f.write('l'+str(ti)+': line = ('+', '.join(line)+');\n')
        linelist.append('l'+str(ti))
        ti += 1
    # need to define the period before making sampler planes
    f.write('lattice: line = ('+', '.join(linelist)+');\n')
    f.write('use, period=lattice;\n')
    f.close()

    #write samplers
    # if less than 10 samplers, just put in main file
    if len(machine.samplers) > 10:
        f = open(fn_samplers,'w')
        files.append(fn_samplers)
        f.write(timestring)
        f.write('! pybdsim.Builder Lattice \n')
        f.write('! SAMPLER DEFINITION\n\n')
        for sampler in machine.samplers:
            f.write(str(sampler))
        f.close()

    # write beam 
    f = open(fn_beam,'w') 
    files.append(fn_beam)
    f.write(timestring) 
    f.write('! pybdsim.Builder \n')
    f.write('! BEAM DEFINITION \n\n')
    f.write(machine.beam.ReturnBeamString())
    f.close()

    # write options - only if specified
    if machine.options != None:
        f = open(fn_options,'w')
        files.append(fn_options)
        f.write(timestring) 
        f.write('! pybdsim.Builder \n')
        f.write('! OPTIONS DEFINITION \n\n')
        f.write(machine.options.ReturnOptionsString())
        f.close()

    # write main file
    f = open(fn_main,'w')
    f.write(timestring)
    f.write('! pybdsim.Builder Lattice \n')
    f.write('! number of elements = ' + str(len(machine.elements)) + '\n')
    f.write('! total length       = ' + str(machine.length) + ' m\n\n')
    
    for fn in files:
        fn = fn.split('/')[-1]
        f.write('include '+fn+';\n')
    if len(machine.samplers) <= 10:
        for sampler in machine.samplers:
            f.write(str(sampler))
    f.close()

    #user feedback
    print('Lattice written to:')
    for fn in files:
        print(fn)
    print('All included in main file: \n',fn_main)


def GenerateSamplersFromBDSIMSurvey(surveyfile,outputfilename,excludesamplers=True):
    """
    Create a gmad file with samplers for all the elements in a beamline
    as described by the survey outline from bdsim
    
    bdsim --file=mylattice.gmad --outline=survey.dat --outline_type=survey

    excludesamplers - bool - exclude any existing samplers

    """
    a = Data.Load(surveyfile)
    samplers = []
    for name in a.Name():
        if ('ampler' in name) and excludesamplers:
            pass
        else:
            samplers.append(Sampler(name))

    #write the output
    f = open(outputfilename,'w')
    for sampler in samplers:
        f.write(sampler.__repr__())
    f.close()
