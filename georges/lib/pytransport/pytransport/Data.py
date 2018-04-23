# pytransport.Data - Data loaders and data containers.
# Version 1.0
# W. Shields and J. Snuverink
# william.shields.2010@live.rhul.ac.uk

"""
Data

Data containers used in converting from Transport to gmad or madx.

Classes:
BDSData - a list of data read from Transport files.
ConversionData - a class for holding data during conversion.

"""

import numpy as _np
import os as _os
from scipy import constants as _con
import copy

_useRootNumpy = True

try:
    import root_numpy as _rnp
except ImportError:
    _useRootNumpy = False
    pass


def _Load(filepath):
    extension = filepath.split('.')[-1]
    if not _os.path.isfile(filepath):
        raise IOError("File does not exist")
    if ("elosshist" in filepath) or (".hist" in filepath):
        return _LoadAsciiHistogram(filepath)
    elif extension == 'root':
        try:
            return _LoadRoot(filepath)
        except NameError:
            # raise error rather than return None, saves later scripting errors.
            raise IOError('Root loader not available.')
    else:
        raise IOError("Unknown file type - not BDSIM data")


def _LoadAsciiHistogram(filepath):
    data = BDSData()
    f = open(filepath, 'r')
    for i, line in enumerate(f):
        # first line is header (0 counting)
        if i == 1:
            names, units = _ParseHeaderLine(line)
            for name, unit in zip(names, units):
                data._AddProperty(name, unit)
        elif "nderflow" in line:
            data.underflow = float(line.strip().split()[1])
        elif "verflow" in line:
            data.overflow = float(line.strip().split()[1])
        elif i >= 4:
            data.append(tuple(map(float, line.split())))
    f.close()
    return data


def _LoadRoot(filepath):
    if not _useRootNumpy:
        raise IOError("root_numpy not available - can't load ROOT file")
    data = BDSData()
    trees = _rnp.list_trees(filepath)

    if 'optics' in trees:
        branches = _rnp.list_branches(filepath, 'optics')
        treedata = _rnp.root2array(filepath, 'optics')
    elif 'orbit' in trees:
        branches = _rnp.list_branches(filepath, 'orbit')
        treedata = _rnp.root2array(filepath, 'orbit')
    else:
        raise IOError("This file doesn't have the required tree 'optics'.")
    for element in range(len(treedata[branches[0]])):
        elementlist = []
        for branch in branches:
            if element == 0:
                data._AddProperty(branch)
            elementlist.append(treedata[branch][element])
        data.append(elementlist)
    return data


def _ParseHeaderLine(line):
    names = []
    units = []
    for word in line.split():
        if word.count('[') > 0:
            names.append(word.split('[')[0])
            units.append(word.split('[')[1].strip(']'))
        else:
            names.append(word)
            units.append('NA')
    return names, units


class BDSData(list):
    """
    General class representing simple 2 column data.

    Inherits python list.  It's a list of tuples with extra columns of 'name' and 'units'.
    """
    def __init__(self, *args, **kwargs):
        list.__init__(self, *args, **kwargs)
        self.units   = []
        self.names   = []
        self.columns = self.names

    def __getitem__(self, index):
        return dict(zip(self.names, list.__getitem__(self, index)))

    def GetItemTuple(self, index):
        """
        Get a specific entry in the data as a tuple of values rather than a dictionary.
        """
        return list.__getitem__(self, index)
        
    def _AddMethod(self, variablename):
        """
        This is used to dynamically add a getter function for a variable name.
        """
        def GetAttribute():
            if self.names.count(variablename) == 0:
                raise KeyError(variablename+" is not a variable in this data")
            ind = self.names.index(variablename)
            return _np.array([event[ind] for event in self])
        setattr(self, variablename, GetAttribute)

    def ConcatenateMachine(self, *args):
        """
        This is used to concatenate machines.
        """
        # Get final position of the machine (different param for survey)
        lastSpos = self.GetColumn('S')[-1]
        
        for machine in args:
            if isinstance(machine, _np.str):
                machine = _Load(machine)
        
            # check names sets are equal
            if len(set(self.names).difference(set(machine.names))) != 0:
                raise AttributeError("Cannot concatenate machine, variable names do not match")

            if self.names.count('S') != 0:
                sind = self.names.index('S')
            else:
                raise KeyError("S is not a variable in this data")
        
            # Have to convert each element to a list as tuples can't be modified
            for index in range(len(machine)):
                element = machine.GetItemTuple(index)
                elementlist = list(element)
                elementlist[sind] += lastSpos
                
                self.append(tuple(elementlist))

            lastSpos += machine.GetColumn('S')[-1]

    def _AddProperty(self, variablename, variableunit='NA'):
        """
        This is used to add a new variable and hence new getter function
        """
        self.names.append(variablename)
        self.units.append(variableunit)
        self._AddMethod(variablename)

    def _DuplicateNamesUnits(self, bdsdata2instance):
        d = bdsdata2instance
        for name, unit in zip(d.names, d.units):
            self._AddProperty(name, unit)

    def MatchValue(self, parametername, matchvalue, tolerance):
        """
        This is used to filter the instance of the class based on matching
        a parameter withing a certain tolerance.

        >>> a = pytransport.Data.Load("myfile.txt")
        >>> a.MatchValue("S",0.3,0.0004)
        
        this will match the "S" variable in instance "a" to the value of 0.3
        within +- 0.0004.

        You can therefore used to match any parameter.

        Return type is BDSAsciiData
        """
        if hasattr(self, parametername):
            a = BDSData()                 #build bdsdata2
            a._DuplicateNamesUnits(self)  #copy names and units
            pindex = a.names.index(parametername)
            filtereddata = [event for event in self if abs(event[pindex] - matchvalue) <= tolerance]
            a.extend(filtereddata)
            return a
        else:
            print "The parameter: ", parametername, " does not exist in this instance"

    def Filter(self, booleanarray):
        """
        Filter the data with a booleanarray.  Where true, will return
        that event in the data.

        Return type is BDSData
        """
        a = BDSData()
        a._DuplicateNamesUnits(self)
        a.extend([event for i, event in enumerate(self) if booleanarray[i]])
        return a

    def NameFromNearestS(self, S):
        i = self.IndexFromNearestS(S)
        if not hasattr(self, "Name"):
            raise ValueError("This file doesn't have the required column Name")
        return self.Name()[i]
    
    def IndexFromNearestS(self, S):
        """
        IndexFromNearestS(S) 

        return the index of the beamline element clostest to S 

        Only works if "SStart" column exists in data
        """
        # check this particular instance has the required columns for this function
        if not hasattr(self, "SStart"):
            raise ValueError("This file doesn't have the required column SStart")
        if not hasattr(self, "Arc_len"):
            raise ValueError("This file doesn't have the required column Arc_len")
        s = self.SStart()
        l = self.Arc_len()

        # iterate over beamline and record element if S is between the
        # sposition of that element and then next one
        # note madx S position is the end of the element by default
        ci = [i for i in range(len(self) - 1) if (S > s[i] and S < s[i]+l[i])]
        try:
            ci = ci[0] # return just the first match - should only be one
        except IndexError:
            # protect against S positions outside range of machine
            if S > s[-1]:
                ci -= 1
            else:
                ci = 0
        # check the absolute distance to each and return the closest one
        # make robust against s positions outside machine range
        return ci

    def GetColumn(self, columnstring):
        """
        Return a numpy array of the values in columnstring in order
        as they appear in the beamline
        """
        if columnstring not in self.columns:
            raise ValueError("Invalid column name")
        ind = self.names.index(columnstring)
        return _np.array([element[ind] for element in self])

    def __repr__(self):
        s = ''
        s += 'pytransport.Data.BDSData instance\n'
        s += str(len(self)) + ' entries'
        return s


class ConversionData:
    """
    Class used as data container object in Transport2Gmad / Transport2Madx conversion.
    Required input:
    - inputfile: string, inputfile name
    - machine: either pybdsim.Builder.Machine or pymadx.Builder.Machine instance.

    Note: if used as a holder for conversion to gmad, options must be supplied a pybdsim.Options.Options instance.

    This class will hold ALL conversion related data, some stored in member variables which are separate containers
    for: conversion related properties (user input arguments), beam properties, and machine properties.
    """
    def __init__(self,
                 inputfile,
                 machine,
                 options       = None,  # None as madx has no options class. gmad needs pybdsim.Options passing in
                 particle      = 'proton',
                 debug         = False,
                 distrType     = 'gauss',
                 gmad          = True,
                 gmadDir       = 'gmad',
                 madx          = False,
                 madxDir       = 'madx',
                 auto          = True,
                 dontSplit     = False,
                 keepName      = False,
                 combineDrifts = False,
                 outlog        = True):

        if particle == 'proton':
            p_mass = _con.proton_mass * (_con.c ** 2 / _con.e) / 1e9  # Particle masses in same unit as TRANSPORT (GeV)
        elif particle == 'e-' or particle == 'e+':
            p_mass = _con.electron_mass * (_con.c ** 2 / _con.e) / 1e9
        else:
            p_mass = 1

        # pytransport data container classes. Split into different container classes so that this class is not
        # overloaded with member variables. Could be stored as a dictionary but given their widespread use in
        # conversion, access is just cleaner.
        self.convprops = _conversionProps(inputfile, particle, debug, gmad, gmadDir, madx, madxDir,
                                                   auto, dontSplit, keepName, combineDrifts, outlog)
        self.beamprops = _beamprops(p_mass)
        self.beamprops.distrType = distrType
        self.machineprops = _machineprops()

        # can't verify correct data type so at least check that it's not None if converting to gmad.
        if (options is None) and gmad:
            raise TypeError("options must be a pybdsim.Options.Options instance when converting to gmad")
        else:
            self.options = options

        # the gmad/madx machine and beam that will be written.
        self.machine = machine
        self.beam = self.machine.beam
        self.beam['offsetSampleMean'] = 0

        # make a copy of the empty machine. Copy needed in case machine is split and a new machine is needed.
        self._machineCopy = copy.deepcopy(self.machine)

        # initialise registries
        self.ElementRegistry = _Registry()
        self.FitRegistry = _Registry()

        self.units = {  # Default TRANSPORT units
            'x': 'cm',
            'xp': 'mrad',
            'y': 'cm',
            'yp': 'mrad',
            'bunch_length': 'cm',
            'momentum_spread': 'pc',
            'element_length': 'm',
            'magnetic_fields': 'kG',
            'p_egain': 'GeV',  # Momentum / energy gain during acceleration.
            'bend_vert_gap': 'cm',  # Vertical half-gap in dipoles
            'pipe_rad': 'cm',
            'beta_func': 'm',
            'emittance': 'mm mrad'
        }
        self.scale = {
            'p': 1e-12,
            'n': 1e-9,
            'u': 1e-6,
            'm': 1e-3,
            'c': 1e-2,
            'k': 1e+3,
            'K': 1e+3,  # Included both cases of k just in case.
            'M': 1e+6,
            'G': 1e+9,
            'T': 1e+12
        }

        self.accstart = []  # An index of the start of acceleration elements.
        self.data = []  # A list that will contain arrays of the element data
        self.filedata = []  # A list that will contain the raw strings from the input file

    def AddOptions(self):
        """
        Function to set the Options for a BDSIM machine.
        """
        self.options.SetPhysicsList(physicslist='em')
        self.options.SetBeamPipeRadius(beampiperadius=self.machineprops.beampiperadius,
                                       unitsstring=self.units['pipe_rad'])
        self.options.SetOuterDiameter(outerdiameter=0.5, unitsstring='m')
        self.options.SetTunnelRadius(tunnelradius=1, unitsstring='m')
        self.options.SetBeamPipeThickness(bpt=5, unitsstring='mm')
        self.options.SetSamplerDiameter(radius=1, unitsstring='m')
        self.options.SetStopSecondaries(stop=True)
        self.options.SetIncludeFringeFields(on=True)

        self.machine.AddOptions(self.options)

    def AddBeam(self):
        """
        Function to prepare the beam and add to the machine.
        """
        # convert energy to GeV (madx only handles GeV)
        energy_in_gev = self.beamprops.tot_energy * self.scale[self.units['p_egain'][0]] / 1e9
        self.beamprops.tot_energy = energy_in_gev

        self.beam.SetParticleType(self.convprops.particle)
        self.beam.SetEnergy(energy=self.beamprops.tot_energy, unitsstring='GeV')

        if self.convprops.gmadoutput:
            # set gmad parameters depending on distribution
            if self.beamprops.distrType == 'gausstwiss':
                self.beam.SetDistributionType(self.beamprops.distrType)
                self.beam.SetBetaX(self.beamprops.betx)
                self.beam.SetBetaY(self.beamprops.bety)
                self.beam.SetAlphaX(self.beamprops.alfx)
                self.beam.SetAlphaY(self.beamprops.alfy)
                self.beam.SetEmittanceX(self.beamprops.emitx, unitsstring='mm')
                self.beam.SetEmittanceY(self.beamprops.emity, unitsstring='mm')
                self.beam.SetSigmaE(self.beamprops.SigmaE)
                self.beam.SetSigmaT(self.beamprops.SigmaT)

            else:
                self.beam.SetDistributionType(self.beamprops.distrType)
                self.beam.SetSigmaX(self.beamprops.SigmaX, unitsstring=self.units['x'])
                self.beam.SetSigmaY(self.beamprops.SigmaY, unitsstring=self.units['y'])
                self.beam.SetSigmaXP(self.beamprops.SigmaXP, unitsstring=self.units['xp'])
                self.beam.SetSigmaYP(self.beamprops.SigmaYP, unitsstring=self.units['yp'])
                self.beam.SetSigmaE(self.beamprops.SigmaE)
                self.beam.SetSigmaT(self.beamprops.SigmaT)

            # set beam offsets in gmad if non zero
            if self.beamprops.X0 != 0:
                self.beam.SetX0(self.beamprops.X0, unitsstring=self.units['x'])
            if self.beamprops.Y0 != 0:
                self.beam.SetY0(self.beamprops.Y0, unitsstring=self.units['y'])
            if self.beamprops.Z0 != 0:
                self.beam.SetZ0(self.beamprops.Z0, unitsstring=self.units['z'])

        elif self.convprops.madxoutput:
            # calculate betas and emittances regardless for madx beam
            try:
                self.beamprops.betx = self.beamprops.SigmaX / self.beamprops.SigmaXP
            except ZeroDivisionError:
                self.beamprops.betx = 0
            try:
                self.beamprops.bety = self.beamprops.SigmaY / self.beamprops.SigmaYP
            except ZeroDivisionError:
                self.beamprops.bety = 0
                self.beamprops.emitx = self.beamprops.SigmaX * self.beamprops.SigmaXP / 1000.0
                self.beamprops.emity = self.beamprops.SigmaY * self.beamprops.SigmaYP / 1000.0

            # set madx beam
            self.beam.SetDistributionType('madx')
            self.beam.SetBetaX(self.beamprops.betx)
            self.beam.SetBetaY(self.beamprops.bety)
            self.beam.SetAlphaX(self.beamprops.alfx)
            self.beam.SetAlphaY(self.beamprops.alfy)
            self.beam.SetEmittanceX(self.beamprops.emitx / 1000)
            self.beam.SetEmittanceY(self.beamprops.emity / 1000)
            self.beam.SetSigmaE(self.beamprops.SigmaE)
            self.beam.SetSigmaT(self.beamprops.SigmaT)

        self.machine.AddBeam(self.beam)

    def _NewMachines(self):
        """
        Delete the machine and set to be the empty machine copied at class instantiation.
        """
        del self.machine
        self.machine = self._machineCopy


class _beamprops:
    """
    A class containing the properties of the beam distribution.
    """
    def __init__(self, p_mass=938.272):
        # beam properties that are updated along the lattice
        self.momentum = 0
        self.k_energy = 0
        self.tot_energy_current = p_mass
        self.gamma = 1
        self.beta = 0
        self.brho = 0
        # beam properties that are from the initial beam and fixed
        self.mass = p_mass
        self.tot_energy = p_mass  # initial energy
        self.SigmaX = 0
        self.SigmaY = 0
        self.SigmaXP = 0
        self.SigmaYP = 0
        self.SigmaE = 0
        self.SigmaT = 0
        self.X0 = 0
        self.Y0 = 0
        self.Z0 = 0
        self.T0 = 0
        self.Xp0 = 0
        self.Yp0 = 0
        self.betx = 0
        self.alfx = 0
        self.bety = 0
        self.alfy = 0
        self.dx = 0
        self.dy = 0
        self.emitx = 0
        self.emity = 0
        self.distrType = 'gauss'


class _machineprops:
    """
    A class containing the number of elements and angular properties (i.e bending direction)
    """
    def __init__(self):
        self.benddef        = True  # True = dipole defined by 4. L B n. False = dipole defined by 4. L angle n.
        self.bending        = 1     # +VE = bends to the right for positive particles
        self.angle          = 0     # dipole rotation angle
        self.drifts         = 0     # nr of drifts
        self.dipoles        = 0
        self.rf             = 0
        self.quads          = 0
        self.sextus         = 0
        self.transforms     = 0
        self.solenoids      = 0
        self.collimators    = 0
        self.beampiperadius = 20
        self.fringeIntegral = 0  # global value for all subsequent fringe fields until set otherwise
        self.dipoleVertAper = 0
        self.apertureType   = 'circular'
        self._totalAccVoltage = 0
        self._e_gain_prev   = 0


class _conversionProps:
    """
    A class containing the settings for the conversion.
    """
    def __init__(self, inputfile,
                 particle      = 'proton',
                 debug         = False,
                 gmad          = True,
                 gmadDir       = 'gmad',
                 madx          = False,
                 madxDir       = 'madx',
                 auto          = True,
                 dontSplit     = False,
                 keepName      = False,
                 combineDrifts = False,
                 outlog        = True):

        self.debug = debug
        self.outlog = outlog
        self.typeCode6IsTransUpdate = True  # Definition of type code 6, true is transform update, false is collimator
        self.isAccSequence = False  # Definition of type code 11 is not an accelerator sequence

        # Automatic writing and machine splitting
        self.auto = auto
        self.dontSplit = dontSplit

        # beam definition
        self.particle = particle
        self.beamdefined = False
        self.correctedbeamdef = False

        # File input and output
        self.file = inputfile
        self.fileloaded = False
        self.gmadoutput = gmad
        self.gmadDir = gmadDir
        self.madxoutput = madx
        self.madxDir = madxDir
        self.numberparts = -1

        # transport optics output is modified to be a single line
        self.singleLineOptics = False
        self.keepName = keepName
        self.combineDrifts = combineDrifts


class _Registry:
    """
    A class used as a registry during conversion.
    """
    def __init__(self):
        self.elements = []
        self.names = []
        self.lines = []
        self.length = []
        self._uniquenames = []
        self._totalLength = 0

    def AddToRegistry(self, linedict, line):
        if not isinstance(linedict, dict):
            raise TypeError("Added element is not a Dictionary")
        self.elements.append(linedict)
        self.names.append(linedict['name'])
        if not linedict['name'] in self._uniquenames:
            self._uniquenames.append(linedict['name'])

        self.lines.append(line)
        # Cumulative length
        length = round(linedict['length'], 5)
        if len(self.length) > 0:
            self.length.append(length + self._totalLength)
        else:
            self.length.append(length)
        self._totalLength += length

    def GetElementIndex(self, name):
        elenums = []
        if name not in self.names:
            return elenums
        else:
            # Add all elements of the same name as a single element may be
            # used multiple times.
            for index, elename in enumerate(self.names):
                if elename == name:
                    elenums.append(index)
            return elenums

    def GetElement(self, name):
        elenum = self.GetElementIndex(name)
        if isinstance(elenum, list):
            elementList = []
            for num in elenum:
                elementList.append(self.elements[num])
            return elementList
        else:
            return self.elements[elenum]

    def GetElementEndSPosition(self, name):
        elenum = self.GetElementIndex(name)
        if isinstance(elenum, list):
            elementList = []
            for num in elenum:
                elementList.append(self.length[num])
            return elementList
        else:
            return self.length[elenum]

    def GetElementStartSPosition(self, name):
        elenum = self.GetElementIndex(name)
        endS = self.GetElementEndSPosition(name)

        if isinstance(elenum, list):
            elementList = []
            for index, num in enumerate(elenum):
                element = self.elements[num]
                length = element['length']
                startS = endS[index] - length
                elementList.append(round(startS, 5))
            return elementList
        else:
            element = self.elements[elenum]
            length = element['length']
            startS = endS - length
            return round(startS, 5)

    def UpdateLength(self, linedict):
        """
        Function to increases the machines length, but does not add element data.
        This is so the S positions of named elements in the fitting registry can
        be calculated correctly.
        """
        if not isinstance(linedict, dict):
            raise TypeError("Added element is not a Dictionary")
        self._totalLength += linedict['length']
