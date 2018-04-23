# pybdsim.Data - output data loader for pybdsim
# Version 1.0
# L. Nevay, S.Boogert
# laurie.nevay@rhul.ac.uk

"""
Output

Read bdsim output

Classes:
Data - read various output files


"""
import numpy as _np
from . import Constants as _Constants
from . import _General
import os as _os

useRootNumpy = True

try:
    import root_numpy as _rnp
except ImportError:
    useRootNumpy = False
    pass

def Load(filepath):
    extension = filepath.split('.')[-1]
    if not _os.path.isfile(filepath):
        raise IOError("File does not exist")
    if ("elosshist" in filepath) or (".hist" in filepath):
        return _LoadAsciiHistogram(filepath)
    elif "eloss" in filepath:
        return _LoadAscii(filepath)
    elif extension == 'txt':
        return _LoadAscii(filepath)
    elif extension == 'root':
        try:
            return _LoadRoot(filepath)
        except NameError:
            #raise error rather than return None, saves later scripting errors.
            raise IOError('Root loader not available.')
    elif extension == 'dat':
        print('.dat file - trying general loader')
        try:
            return _LoadAscii(filepath)
        except:
            print("Didn't work")
            raise IOError("Unknown file type - not BDSIM data")
    else:
        raise IOError("Unknown file type - not BDSIM data")

def _LoadAscii(filepath):
    data = BDSAsciiData()
    f = open(filepath, 'r')
    for i, line in enumerate(f):
        if line.startswith("#"):
            pass
        elif i == 1:
        # first line is header
            names,units = _ParseHeaderLine(line)
            for name,unit in zip(names,units):
                data._AddProperty(name,unit)
        else:
            #this tries to cast to float, but if not leaves as string
            data.append(tuple(map(_General.Cast,line.split())))
    f.close()
    return data

def _LoadAsciiHistogram(filepath):
    data = BDSAsciiData()
    f = open(filepath,'r')
    for i, line in enumerate(f):
        # first line is header (0 counting)
        if i == 1:
            names,units = _ParseHeaderLine(line)
            for name,unit in zip(names,units):
                data._AddProperty(name,unit)
        elif "nderflow" in line:
            data.underflow = float(line.strip().split()[1])
        elif "verflow" in line:
            data.overflow  = float(line.strip().split()[1])
        elif i >= 4:
            data.append(tuple(map(float,line.split())))
    f.close()
    return data

def _LoadRoot(filepath):
    if not useRootNumpy:
        raise IOError("root_numpy not available - can't load ROOT file")
    data = BDSAsciiData()
    trees = _rnp.list_trees(filepath)

    if 'optics' in trees:
        branches = _rnp.list_branches(filepath,'optics')
        treedata = _rnp.root2array(filepath,'optics')
    elif 'orbit' in trees:
        branches = _rnp.list_branches(filepath, 'orbit')
        treedata = _rnp.root2array(filepath, 'orbit')
    else:
        raise IOError("This file doesn't have the required tree 'optics'.")
    for element in range(len(treedata[branches[0]])):
        elementlist=[]
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
                

class BDSAsciiData(list):
    """
    General class representing simple 2 column data.

    Inherits python list.  It's a list of tuples with extra columns of 'name' and 'units'.
    """
    def __init__(self, *args, **kwargs):
        list.__init__(self, *args, **kwargs)
        self.units   = []
        self.names   = []
        self.columns = self.names

    def __getitem__(self,index):
        return dict(zip(self.names,list.__getitem__(self,index)))

    def GetItemTuple(self,index):
        """
        Get a specific entry in the data as a tuple of values rather than a dictionary.
        """
        return list.__getitem__(self,index)
        
    def _AddMethod(self, variablename):
        """
        This is used to dynamically add a getter function for a variable name.
        """
        def GetAttribute():
            if self.names.count(variablename) == 0:
                raise KeyError(variablename+" is not a variable in this data")
            ind = self.names.index(variablename)
            return _np.array([event[ind] for event in self])
        setattr(self,variablename,GetAttribute)

    def ConcatenateMachine(self,*args):
        """
        This is used to concatenate machines.
        """
        #Get final position of the machine (different param for survey)
        if _General.IsSurvey(self):
            lastSpos = self.GetColumn('SEnd')[-1]
        else:
            lastSpos = self.GetColumn('S')[-1]
        
        for machine in args:
            if isinstance(machine,_np.str):
                machine = Load(machine)
        
            #check names sets are equal
            if len(set(self.names).difference(set(machine.names))) != 0:
                raise AttributeError("Cannot concatenate machine, variable names do not match")
        
            #surveys have multiple s positions per element
            if _General.IsSurvey(machine):
                sstartind = self.names.index('SStart')
                smidind = self.names.index('SMid')
                sendind = self.names.index('SEnd')
            elif self.names.count('S') != 0:
                sind = self.names.index('S')
            else:
                raise KeyError("S is not a variable in this data")
        
            #Have to convert each element to a list as tuples can't be modified
            for index in range(len(machine)):
                element = machine.GetItemTuple(index)
                elementlist = list(element)
                
                #update the elements S position
                if _General.IsSurvey(machine):
                    elementlist[sstartind] += lastSpos
                    elementlist[smidind] += lastSpos
                    elementlist[sendind] += lastSpos
                else:
                    elementlist[sind] += lastSpos
                
                self.append(tuple(elementlist))
                
            #update the total S position.
            if _General.IsSurvey(machine):
                lastSpos += machine.GetColumn('SEnd')[-1]
            else:
                lastSpos += machine.GetColumn('S')[-1]


    def _AddProperty(self,variablename,variableunit='NA'):
        """
        This is used to add a new variable and hence new getter function
        """
        self.names.append(variablename)
        self.units.append(variableunit)
        self._AddMethod(variablename)

    def _DuplicateNamesUnits(self,bdsasciidata2instance):
        d = bdsasciidata2instance
        for name,unit in zip(d.names,d.units):
            self._AddProperty(name,unit)

    def MatchValue(self,parametername,matchvalue,tolerance):
        """
        This is used to filter the instance of the class based on matching
        a parameter withing a certain tolerance.

        >>> a = pybdsim.Data.Load("myfile.txt")
        >>> a.MatchValue("S",0.3,0.0004)
        
        this will match the "S" variable in instance "a" to the value of 0.3
        within +- 0.0004.

        You can therefore used to match any parameter.

        Return type is BDSAsciiData
        """
        if hasattr(self,parametername):
            a = BDSAsciiData()            #build bdsasciidata2
            a._DuplicateNamesUnits(self)   #copy names and units
            pindex = a.names.index(parametername)
            filtereddata = [event for event in self if abs(event[pindex]-matchvalue)<=tolerance]
            a.extend(filtereddata)
            return a
        else:
            print("The parameter: ",parametername," does not exist in this instance")

    def Filter(self,booleanarray):
        """
        Filter the data with a booleanarray.  Where true, will return
        that event in the data.

        Return type is BDSAsciiData
        """
        a = BDSAsciiData()
        a._DuplicateNamesUnits(self)
        a.extend([event for i,event in enumerate(self) if booleanarray[i]])
        return a

    def NameFromNearestS(self,S):
        i = self.IndexFromNearestS(S)
        if not hasattr(self,"Name"):
            raise ValueError("This file doesn't have the required column Name")
        return self.Name()[i]
    
    def IndexFromNearestS(self,S) : 
        """
        IndexFromNearestS(S) 

        return the index of the beamline element clostest to S 

        Only works if "SStart" column exists in data
        """
        #check this particular instance has the required columns for this function
        if not hasattr(self,"SStart"):
            raise ValueError("This file doesn't have the required column SStart")
        if not hasattr(self,"Arc_len"):
            raise ValueError("This file doesn't have the required column Arc_len")
        s = self.SStart()
        l = self.Arc_len()

        #iterate over beamline and record element if S is between the
        #sposition of that element and then next one
        #note madx S position is the end of the element by default
        ci = [i for i in range(len(self)-1) if (S > s[i] and S < s[i]+l[i])]
        try:
            ci = ci[0] #return just the first match - should only be one
        except IndexError:
            #protect against S positions outside range of machine
            if S > s[-1]:
                ci =-1
            else:
                ci = 0
        #check the absolute distance to each and return the closest one
        #make robust against s positions outside machine range
        return ci

    def GetColumn(self,columnstring):
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
        s += 'pybdsim.Data.BDSAsciiData instance\n'
        s += str(len(self)) + ' entries'
        return s
