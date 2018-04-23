# pybdsim.Analysis - analysis scripts for bdsim output
# Version 1.0
# L. Nevay, S.T.Boogert
# laurie.nevay@rhul.ac.uk

#this module does not feature directly in the pybdsim import
#but instead acts as a holder for the analysis class to
#keep the __init__ file very clean.

import Data
import numpy as _np
import Plot
import Constants

class Analysis:
    """
    Analysis class for bdsim output

    one instance of the class per output file

    Analysis('../../path/to/my/output.txt')
    
    It will also look for '../../path/to/my/output.eloss.txt' 
    beside the output file

    """
    def __init__(self,filepath):
        self.data,self.dataarray,self.keyslist = Data.LoadOld(filepath)
        self.filepath = filepath
        self.filename = filepath.split('/')[-1]
        self.plots     = []
        self._plottitle = self.filename.split('.')[0]
        self.ndata     = len(self.data['X'])
        
    def GroupBy(self,variable='Z'):
        """
        GroupBy(variable='Z')
        
        create instance.datagrouped dictionary
        finds unique values of variable and groups all the 
        data that has that value of variable

        e.g.
        GroupBy()   # default is 'Z'
        instance.datagrouped is dict with:
        0.000: array of (nparticles x allother variables)
        1.202: similar but different number of particles maybe
        ... etc
        """
        self.datagrouped = Data.AsciiData()
        
        #find unique values of variable
        uniquevalues = sorted(list(set(_np.round(self.data[variable],2))))

        #remove the variable from the subset
        #find it's index in keys list
        indexofvariabletoremove = self.keyslist.index(variable)
        for value in uniquevalues:
            mask      = _np.round(self.data[variable],2) == value
            dcopy     = self.dataarray[mask]
            dcopydict = dict(zip(self.keyslist,[dcopy[:,i] for i in range(_np.shape(dcopy)[1])]))

            #dcopydict = Data.AsciiData(zip(self.data.keyslist,[dcopy[:,i] for i in range(_np.shape(dcopy)[1])]))
            dcopydict['nparticles'] = _np.shape(dcopy)[0]
            self.datagrouped[value] = dcopydict
        self.keysgrouped = list(_np.sort(self.datagrouped.keys()))

    def GenerateSigmas(self):
        if hasattr(self,'datagrouped') == False:
            self.GroupBy()
        
        z = self.keysgrouped
        sx,sy = [],[]

        for i in range(len(self.keysgrouped)):
            sx.append(_np.std(self.datagrouped[self.keysgrouped[i]]['X']))
            sy.append(_np.std(self.datagrouped[self.keysgrouped[i]]['Y']))
        self.simpledata = {'sx':sx,'sy':sy,'z':z}

    def SortBy(self,variable='Z'):
        pass

    def PlotXXp(self):
        p1 = Plot.PlotXY(self.data['X'],self.data['Xp'],'X ('+self.data.units['X']+')','Xp ('+self.data.units['Xp']+')',self._plottitle)
        self.plots.append(p1)

    def PlotYYp(self):
        p1 = Plot.PlotXY(self.data['Y'],self.data['Yp'],'Y ('+self.data.units['Y']+')','Yp ('+self.data.units['Y']+')',self._plottitle)
        self.plots.append(p1)

    def PlotSingleVariable(self,variable,xlabel='',ylabel=''):
        p1 = Plot.PlotY(self.data[variable],xlabel,ylabel,self._plottitle)
        self.plots.append(p1)

    def ControlPlotAll(self):
        p1 = Plot.ControlPlot(self.data,self._plottitle)
        self.plots.append(p1)
        
    def PlotEnergyHist(self,nbins=20,**kwargs):
        p1 = Plot.Histogram(self.data['E'],'Energy ('+self.data.units['E']+')',self._plottitle,nbins,**kwargs)
        self.plots.append(p1)

    

