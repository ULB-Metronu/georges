import ROOT as _ROOT

class Run :
    def __init__(self, filename):
        self._filename = filename
        self._rootFile = _ROOT.TFile(filename)
        self._tree     = self._rootFile.Get("Run")
        self._entryInfo = _ROOT.BDSOutputROOTEventRunInfo()
        self._entryHist = _ROOT.BDSOutputROOTEventHistograms()
        self._tree.SetBranchAddress("Info.",self._entryInfo)
        self._tree.SetBranchAddress("Histos.",self._entryHist)
        self._tree.GetEntry(0)

