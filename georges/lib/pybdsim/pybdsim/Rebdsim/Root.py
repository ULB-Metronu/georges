"""
Converts ROOT histogram to matplotlib figure
"""

import matplotlib.pyplot as _plt
import matplotlib.cm as _cm
import numpy as _np
import matplotlib.ticker as _ticker
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable



'''
Wrapper classes for ROOT data to be plotted properly in matplotlib
'''

class TTree :
    def __init__(self, r2n_tree):
        '''
        :param r2n_tree: root_numpy converted tree
        '''
        self.data = r2n_tree

    def Plot(self, keyX, keyY, keyYErr = '', opt = '', label = ''):
        if keyYErr == '' :
            _plt.plot(self.data[keyX], self.data[keyY], opt, label = label)
        else :
            _plt.errorbar(self.data[keyX], self.data[keyY], self.data[keyYErr], opt, label = label)

    def Keys(self):
        return self.data.dtype.names

class TH1 :
    # TODO : TH1 Deal with under/overflows properly
    def __init__(self, hist):
        self.hist   = hist

        # extract meta data
        self.name   = hist.GetName()
        self.title  = hist.GetTitle()
        self.labelX = hist.GetXaxis().GetTitle()
        self.labelY = hist.GetYaxis().GetTitle()

        #Extract statistics box data
        self.mean  = hist.GetMean()
        self.stdev = hist.GetStdDev()

        # extract data
        nbinsx   = hist.GetNbinsX()
        self.entries   = hist.GetEntries()
        self.widths   = _np.zeros(nbinsx)
        self.centres  = _np.zeros(nbinsx)
        self.lowedge  = _np.zeros(nbinsx)
        self.highedge = _np.zeros(nbinsx)
        self.contents = _np.zeros(nbinsx)
        self.errors   = _np.zeros(nbinsx)

        for i in range(nbinsx):
            self.widths[i]   = hist.GetXaxis().GetBinWidth(i)
            self.lowedge[i]  = hist.GetBinLowEdge(i+1)
            self.highedge[i] = hist.GetBinLowEdge(i+1)
            self.centres[i]  = hist.GetBinCenter(i+1)
            self.contents[i] = hist.GetBinContent(i+1)
            self.errors[i]   = hist.GetBinError(i+1)

    def Plot(self,opt ='hist', logx=False, logy=False, outputfilename=None, axlabels=None, color='b'):
        '''
        ROOT like plotting options for convenience
        :param opt: hist|e1
        :return:
        '''
        _plt.figure(num=1, figsize=(11,7), facecolor="w", edgecolor="w")

        if opt.find('hist') != -1 :
            self.PlotHist()
        elif opt.find('line') != -1 :
            self.PlotPlot()
        elif opt.find('e1') != -1:
            self.PlotErrorbar(color)

        self._SetLabels(axlabels)
        self._StatBox()

        if logx or logy:
            _LogAxes(logx, logy)

        if outputfilename != None:
            _SaveFigure(outputfilename)

    def PlotPlot(self):
        _plt.plot(self.centres,self.contents)

    def PlotErrorbar(self, edgecolor='none', color='b', label='', colorin='b'):
        p = _plt.plot(self.centres+self.widths/2, self.contents, ls="steps", color=colorin)
        _plt.errorbar(self.centres,self.contents, self.errors, fmt="", color=colorin, label=label)

    def PlotBar(self, edgecolor='none', color='b', label=''):
        _plt.bar(self.lowedge, self.contents, self.widths, edgecolor=edgecolor, color=color, label=label)

    def PlotHist(self, edgecolor='none', color='b', label=''):
        _plt.hist(self.centres, self.lowedge, weights=self.contents, edgecolor=edgecolor, color=color, label=label)

    def _SetLabels(self, axlabels):
        if axlabels:
            _plt.xlabel(axlabels[0])
            _plt.ylabel(axlabels[1])
        else:
            _plt.xlabel(self.labelX)
            _plt.ylabel(self.labelY)

    def _StatBox(self):
        _plt.plot([],[], linestyle="None", label=self.title)
        _plt.plot([],[], linestyle="None", label=r"Entries".ljust(15)+_LegNum(self.entries))
        _plt.plot([],[], linestyle="None", label=r"Mean".ljust(15)+_LegNum(self.mean))
        _plt.plot([],[], linestyle="None", label=r"Std Dev".ljust(14)+_LegNum(self.stdev))
        leg = _plt.legend(loc=1, frameon=True, edgecolor="k",
                    handlelength=0.05, fontsize="xx-small")
        if leg:
            leg.draggable(True,update='bbox')



class TH2 :
    # TODO : TH2 Deal with under/overflows properly
    def __init__(self,hist):
        self.hist   = hist

        # extract meta data
        self.name   = hist.GetName()
        self.title  = hist.GetTitle()
        self.labelX = hist.GetXaxis().GetTitle()
        self.labelY = hist.GetYaxis().GetTitle()

        #Extract statistics box data
        self.meanx  = hist.GetMean(1)
        self.stdevx = hist.GetStdDev(1)
        self.meany  = hist.GetMean(2)
        self.stdevy = hist.GetStdDev(2)

        # extract data
        nbinsx   = hist.GetNbinsX()
        nbinsy   = hist.GetNbinsY()
        self.entries   = hist.GetEntries()
        self.xwidths   = _np.zeros(nbinsx)
        self.xcentres  = _np.zeros(nbinsx)
        self.xlowedge  = _np.zeros(nbinsx)
        self.xhighedge = _np.zeros(nbinsx)
        self.ywidths   = _np.zeros(nbinsy)
        self.ycentres  = _np.zeros(nbinsy)
        self.ylowedge  = _np.zeros(nbinsy)
        self.yhighedge = _np.zeros(nbinsy)
        self.contents = _np.zeros((nbinsx,nbinsy))
        self.errors   = _np.zeros((nbinsx,nbinsy))

        for i in range(nbinsx):
            self.xwidths[i]   = hist.GetXaxis().GetBinWidth(i+1)
            self.xlowedge[i]  = hist.GetXaxis().GetBinLowEdge(i+1)
            self.xhighedge[i] = hist.GetXaxis().GetBinLowEdge(i+2)
            self.xcentres[i]  = hist.GetXaxis().GetBinCenter(i+1)

        for i in range(nbinsy):
            self.ywidths[i]   = hist.GetYaxis().GetBinWidth(i+1)
            self.ylowedge[i]  = hist.GetYaxis().GetBinLowEdge(i+1)
            self.yhighedge[i] = hist.GetYaxis().GetBinLowEdge(i+2)
            self.ycentres[i]  = hist.GetYaxis().GetBinCenter(i+1)

        for i in range(nbinsx) :
            for j in range(nbinsy) :
                self.contents[i,j] = hist.GetBinContent(i+1,j+1)
                self.errors[i,j]   = hist.GetBinError(i+1,j+1)

    def Plot(self, logx=False, logy=False, logz=False, cmap="jet", axlabels=None, outputfilename=None):
        _plt.figure(num=1, figsize=(12,6), facecolor="w", edgecolor="w")

        self.PlotColz(logx,logy, logz, cmap, axlabels) #In the future there may be other plotting options

        #self._SetLabels(axlabels) #TODO: Make the label setting more robust
        self._StatBox()

        if outputfilename != None:
            _SaveFigure(outputfilename)

    def PlotColz(self, logx=False, logy=False, logz=False, cmap="jet", axlabels=None):
        _plt.figure(num=1, figsize=(12,6), facecolor="w", edgecolor="w")
        ax = _plt.gca()
        if axlabels:
            _plt.xlabel(axlabels[0])
            _plt.ylabel(axlabels[1])

        xx, yy = _np.meshgrid(self.xcentres,self.ycentres)
        ccmap = _cm.get_cmap(cmap, 20)
        cts = self.contents #shortcut
        if logz:
            _plt.pcolormesh(xx,yy,cts, cmap=ccmap, norm=LogNorm(vmin=1, vmax=cts.max()))
            pdx = 0.005*(xx.max()-xx.min())
            pdy = 0.005*(yy.max()-yy.min())
            _plt.axis([xx.min()-pdx, xx.max()+pdx, yy.min()-pdy, yy.max()+pdy])
            ax=_plt.gca()
            if logx:
                ax.set_xscale("log")
            if logy:
                ax.set_yscale("log")
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cb = _plt.colorbar(cax=cax,format=_ticker.FuncFormatter(_fmtCbar),ticks=_ticker.LogLocator())
        else:
            _plt.pcolormesh(xx,yy,cts)
            cb = _plt.colorbar(format=_ticker.FuncFormatter(_fmtCbar),ticks=_ticker.LogLocator())

        cb.ax.tick_params(labelsize="xx-small")

    def _SetLabels(self, axlabels):
        if axlabels:
            _plt.xlabel(axlabels[0])
            _plt.ylabel(axlabels[1])
        else:
            _plt.xlabel(self.labelX)
            _plt.ylabel(self.labelY)

    def _StatBox(self):
        _plt.plot([],[], linestyle="None", label=self.title)
        _plt.plot([],[], linestyle="None", label=r"Entries".ljust(15)+_LegNum(self.entries))
        _plt.plot([],[], linestyle="None", label=r"Mean x".ljust(15)+_LegNum(self.meanx))
        _plt.plot([],[], linestyle="None", label=r"Mean y".ljust(15)+_LegNum(self.meany))
        _plt.plot([],[], linestyle="None", label=r"Std Dev x".ljust(14)+_LegNum(self.stdevx))
        _plt.plot([],[], linestyle="None", label=r"Std Dev y".ljust(14)+_LegNum(self.stdevy))
        leg = _plt.legend(bbox_to_anchor=(0.68, 0.62), bbox_transform=_plt.gcf().transFigure,
                          loc="lower left", frameon=True, edgecolor="k", handlelength=0.05,
                          fontsize="xx-small", framealpha=1)
        #if leg:
        #    leg.draggable(True, update='loc')
        #leg.set_zorder(20)

#Unbound utility functions
def _LogAxes(x=False, y=False):
    ax = _plt.gca()
    if x:
        ax.set_xscale("log", nonposy='clip')
    if y:
        ax.set_yscale("log", nonposy='clip')

def _LegNum(number): #TODO: make more elegant with string formatting
    n = number #shortcut
    l = str(n)

    if "e" in l:              #Numbers that are already in scientific notation
        ln  = l.split("e")[0]
        le  = l.split("e")[1]
        try:
            lnr = ln[:7]
        except:
            lnr = ln

        nas = lnr+r"x10$^{"+le+r"}$"

    elif n > 1.e6:             #Large numbers (eg. 2375)
        if l[:-1][-1]==("."):  #Take into account trailing .0s (eg. 2375.0)
            ll = len(l)-2
        else:
            ll = len(l)        #ll is the power of 10 (e.g 3)
        lf = l[0]              #The first number will be the integer base (eg. 2)
        try:
            ls = l[1:6]        #Get the fractional base (e.g 375)
        except:
            ls = l[1:]
        if not ls:             #Or add a zero if no fractional base
            ls = "0"
        nas = lf+r"."+ls+r"x10$^{"+str(ll-1)+r"}$" #Put a decimal point and a power of 10

    elif n < 1.e-6:               #Small numbers (eg. 0.0000456)
        lsp  = l.split(".")[1]    #Get digits after the decimal point
        lr   = lsp.split("0")[-1]         #Get the part with no leading zeros (e.g 456)
        nz   = len(lsp.split("0")[-1])-1  #Count the leading zeros
        lf   = lr[0]                      #Get the desired integer base

        try:
            ls   = lr[1:6]     #Geth the fractional base
        except:
            ls   = lr[1:]
        if not ls:
            ls = "0"           #If no fractional base add a zero
        nas  = lf+r"."+ls+r"x10$^{-"+str(nz+2)+r"}$" #Put decimal point and neg power of 10

    else:
        try:
            nas = l[:7]        #If number is small/big enough print without alteration
        except:
            nas = l

    return nas.ljust(11)


def _fmtCbar(x, pos): #Format in scientific notation and make vals < 1 = 0
    if float(x) == 1.0:
        fst = r"$1$" #For a histogram valuesa smalled that 1 are set to 0
    else:            #Such values are set as dummies to allow log plots
        a, b = '{:.0e}'.format(x).split('e')
        b = int(b)
        fst = r'$10^{{{}}}$'.format(b)
    return fst

def _SaveFigure(outputfilename):
    if '.' in outputfilename:
        outputfilename = outputfilename.split('.')[0]
    _plt.savefig(outputfilename+'.pdf')
    #_plt.savefig(outputfilename+'.png')

