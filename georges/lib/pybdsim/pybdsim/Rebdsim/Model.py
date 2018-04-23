import ROOT as _ROOT

class Model :
    def __init__(self, filename):
        self._filename = filename
        self._rootFile = _ROOT.TFile(filename)
        self._tree     = self._rootFile.Get("Model")
        self._entry = _ROOT.BDSOutputROOTEventModel()
        self._tree.SetBranchAddress("Model.",self._entry)
        self._tree.GetEntry(0)

