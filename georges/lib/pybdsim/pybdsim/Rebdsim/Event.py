import ROOT as _ROOT
import ROOT.gDirectory as _gDirectory
import Root as _Root
import numpy as _np
import platform as _platform


class Event :
    '''
    Converter and convenience methods for Bdsim root output files in python
    '''

    def __init__(self, filename):

        # TODO move over to chains of files
        self._filename = filename
        self._rootFile = _ROOT.TFile(filename)
        print self._rootFile
        if not self._rootFile.IsOpen() :
            raise IOError
        self._tree     = self._rootFile.Get("Event")


        # find numbe of events in file
        self.nevent    = self._tree.GetEntries()

        # set branch addresses etc
        # TODO : Add samplers
        self.info            = _ROOT.BDSOutputROOTEventInfo()
        self.primary         = _ROOT.BDSOutputROOTEventSampler("float")()
        self.eloss           = _ROOT.BDSOutputROOTEventLoss()
        self.primaryFirstHit = _ROOT.BDSOutputROOTEventLoss()
        self.primaryLastHit  = _ROOT.BDSOutputROOTEventLoss()
        self.tunnelHit       = _ROOT.BDSOutputROOTEventLoss()
        self.trajectories    = _ROOT.BDSOutputROOTEventTrajectory()
        self.histos          = _ROOT.BDSOutputROOTEventHistograms()

        self._tree.SetBranchAddress("Info.",self.info)
        self._tree.SetBranchAddress("Primary.",self.primary)
        self._tree.SetBranchAddress("PrimaryFirstHit.",self.primaryFirstHit)
        self._tree.SetBranchAddress("PrimaryLastHit.",self.primaryLastHit)
        self._tree.SetBranchAddress("Eloss.",self.eloss)
        self._tree.SetBranchAddress("Trajectory.",self.trajectories)
        self._tree.SetBranchAddress("Histos.",self.histos)

        # dictionary of samplers
        self.samplerDict = {}

    def enableSampler(self,samplerName):
        self.samplerDict[samplerName] = _ROOT.BDSOutputROOTEventSampler("float")()
        self._tree.SetBranchAddress(samplerName+'.',self.samplerDict[samplerName])

    def getEvent(self,ientry):
        if ientry > -1 and ientry < self.nevent :
            self._tree.GetEntry(ientry)
        else :
            raise IndexError

        '''
        Extract complete event ientry
        :param ientry: name of event branch/leaf
        :return: event object
        '''
        pass

    def getNumpyBranch(self,branchname,selector = ''):
        '''
        Extract array of event tree
        :param branchname: name of event branch/leaf
        :param selector: root section string
        :return: numpy array of data
        '''
        nSelected = self._tree.Draw(branchname,selector,"goff")
        dat       = self._tree.GetV1()
        dat.SetSize(nSelected)
        datArray= _np.array(dat)
        return datArray

    def setEvents(self, selector = '', listname = 'evtList'):
        '''
        :param listname: event list name
        :param selector: selection on the tree
        :return:
        '''
        if selector == '' :
            self._tree.SetEventList(0)

        self._tree.Draw(">>"+listname,selector)
        el = _gDirectory.Get(listname)
        self._tree.SetEventList(el)


    def make1DHistogram(self, command, selector = "", name="hist", title="hist", nbins = -1, xlow=-1., xhigh=1.):
        '''
        Use TTree.Draw to create a histogram, useful when the size of the data cannot be extracted using getNumpyBranch
        :param command:
        :param selector:
        :param name: histogram name
        :param title: histogram title
        :param nbins: number of bins
        :param xlow: histogram lowest edge
        :param xhigh: histogram highest edge
        :return: Root.TH1 object
        '''

        # just in case
        _ROOT.TH1.AddDirectory(True)


        # check for existing histogram
        h = _ROOT.gDirectory.Get(name)
        if h != None :
            _ROOT.gDirectory.Delete(name)

        # make histogram if bins are defined
        if nbins != -1 :
            h = _ROOT.TH1D(name,title,nbins,xlow,xhigh)
        nSelected = self._tree.Draw(command+" >> "+name,selector,"goff")
        rootH = _gDirectory.Get(name)    # root histogram object
        maplH = _Root.TH1(rootH)         # matplotlib histogram object

        return maplH


    def make2DHistogram(self, command, selector = "", name="hist", title="hist",
                        xnbins=-1, xlow=-1., xhigh=1.,
                        ynbins=-1, ylow=-1., yhigh=1.):
        '''
        Use TTree.Draw to create a histogram, useful when the size of the data cannot be extracted using getNumpyBranch
        :param command:
        :param selector:
        :param name: histogram name
        :param title: histogram title
        :param xnbins: x number of bins
        :param xlow: x histogram lowest edge
        :param xhigh: x histogram highest edge
        :param ynbins: y number of bins
        :param ylow: y histogram lowest edge
        :param yhigh: y histogram highest edge
        :return: Root.TH1 object
        '''


        # just in case
        _ROOT.TH2.AddDirectory(True)

        # check for existing histogram
        h = _ROOT.gDirectory.Get(name)
        if h != None :
            _ROOT.gDirectory.Delete(name)

        # make histogram if bins are defined
        if xnbins != -1 and ynbins != -1:
            h = _ROOT.TH2D(name, title, xnbins, xlow, xhigh, ynbins, ylow, yhigh)
        nSelected = self._tree.Draw(command + " >> " + name, selector, "goff")
        rootH = _gDirectory.Get(name)  # root histogram object
        maplH = _Root.TH2(rootH)  # matplotlib histogram object

        return maplH