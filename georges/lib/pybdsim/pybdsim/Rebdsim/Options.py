import ROOT as _ROOT


class Options:
    def __init__(self, filename):
        self._filename = filename
        self._rootFile = _ROOT.TFile(filename)
        self._tree = self._rootFile.Get("Options")
        self._entry = _ROOT.BDSOutputROOTEventOptions()
        self._tree.SetBranchAddress("Options.", self._entry)
        self._tree.GetEntry(0)

