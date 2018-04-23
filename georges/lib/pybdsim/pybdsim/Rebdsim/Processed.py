import ROOT as _ROOT
import root_numpy as _root_numpy
from collections import OrderedDict
import Root as _Root


class Processed:
    '''
    Converter and convenience methods for Rebdsim root output files in python
    '''
    def __init__(self, filename):
        self._filename  = filename
        self._rootFile  = _ROOT.TFile(filename)
        self._rootNumpy = OrderedDict()
        self._rootType  = OrderedDict()
        self._dataDict  = OrderedDict()
        nconvert = self.loadfile();
        if nconvert == 0 :
            print 'Processed> file appears empty'
        else :
            print 'Processed>', nconvert, 'objects converted'

    def loadfile(self):
        '''
        Load and convert root file using root_numpy
        :return: number of converted objects
        '''

        iconverted = 0

        kl = self._rootFile.GetListOfKeys()

        for ki in range(0, kl.GetSize()):
            cname = kl[ki].GetClassName()
            name = kl[ki].GetName()

            if cname == "TDirectoryFile":
                d = self._rootFile.Get(name)
                d.cd()
                # loop over histos in dir
                dkl = d.GetListOfKeys()
                for dki in range(0, dkl.GetSize()):
                    # print "Histogram : ", name + "_" + dkl[dki].GetName()
                    h = d.Get(dkl[dki].GetName())
                    ah = _root_numpy.hist2array(h)
                    self._dataDict[name + "_" + dkl[dki].GetName()] = ah
                    self._rootNumpy[name + "_" + dkl[dki].GetName()] = name + "/" + dkl[dki].GetName()
                    self._rootType[name + "_" + dkl[dki].GetName()] = "TH1"

                    iconverted += 1
            elif cname == "TTree":
                t = self._rootFile.Get(name)
                # print "Tree      : ", name
                at = _root_numpy.tree2array(t)
                self._dataDict[name] = at
                self._rootNumpy[name] = name
                self._rootType[name] = "TTree"
                iconverted += 1

        # self._rootFile.Close()

        return iconverted;

    def keys(self):
        '''

        :return: data keys (ROOT objects) available
        '''
        return self._dataDict.keys()

    def getRootData(self, key):
        '''
        Get raw root data (TH1D, 2D, 3D, TTree)
        :param key: ROOT name of object
        :return: raw root_numpy converted TH1D,2D,3D,TTree
        '''
        return self._rootFile.Get(self._rootNumpy[key])

    def getRootNumpyData(self, key):
        '''
        Get raw root_numpy converted root data (TH1D, 2D, 3D, TTree) using root_numpy
        :param key: ROOT name of object
        :return: raw root_numpy converted TH1D,2D,3D,TTree
        '''

        return self._dataDict[key]

    def getMatplotlibData(self, key):
        '''
        Get root_numpy data converted into python object most useful for plotting with matplotlib
        :param key: ROOT name of object
        :return: root object converted to format most useful for plotting with matplotlib
        '''

        if self._rootType[key] == "TH1" :
            return _Root.TH1(self.getRootData(key))
        elif self._rootType[key] == "TTree" :
            return _Root.TTree(self.getRootNumpyData(key))

