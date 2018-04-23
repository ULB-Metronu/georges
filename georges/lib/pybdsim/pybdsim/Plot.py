#pybdsim plotting tools
# Version 1.0
# L. Nevay, S.T.Boogert
# laurie.nevay@rhul.ac.uk

"""
Useful plots for bdsim output

"""
from . import Data as _Data
import pymadx as _pymadx

import matplotlib as _matplotlib
import matplotlib.pyplot as _plt
import matplotlib.patches as _patches
import numpy as _np
import string as _string

from ._General import CheckItsBDSAsciiData as _CheckItsBDSAsciiData

class _My_Axes(_matplotlib.axes.Axes):
    """
    Inherit matplotlib.axes.Axes but override pan action for mouse.
    Only allow horizontal panning - useful for lattice axes.
    """
    name = "_My_Axes"
    def drag_pan(self, button, key, x, y):
        _matplotlib.axes.Axes.drag_pan(self, button, 'x', x, y) # pretend key=='x'

#register the new class of axes
_matplotlib.projections.register_projection(_My_Axes)

def MadxTfsBetaSimple(tfsfile, title='', outputfilename=None):
    """
    A forward to the pymadx.Plot.PlotTfsBetaSimple function.
    """
    _pymadx.Plot.PlotTfsBetaSimple(tfsfile,title,outputfilename)

def MadxTfsBeta(tfsfile, title='', outputfilename=None):
    """
    A forward to the pymadx.Plot.PlotTfsBeta function.
    """
    _pymadx.Plot.PlotTfsBeta(tfsfile,title,outputfilename)

def AddMachineLatticeToFigure(figure,tfsfile, tightLayout=True):
    """
    A forward to the pymadx.Plot.AddMachineLatticeToFigure function.
    """
    _pymadx.Plot.AddMachineLatticeToFigure(figure, tfsfile, tightLayout)

def ProvideWrappedS(sArray, index):
    s = sArray #shortcut
    smax = s[-1]
    sind = s[index]
    snewa = s[index:]
    snewa = snewa - sind
    snewb = s[:index]
    snewb = snewb + (smax - sind)
    snew  = _np.concatentate((snewa,snewb))
    return snew

def _SetMachineAxesStyle(ax):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

def _PrepareMachineAxes(figure):
    # create new machine axis with proportions 6 : 1
    axmachine = figure.add_subplot(911, projection="_My_Axes")
    _SetMachineAxesStyle(axmachine)
    return axmachine

def _AdjustExistingAxes(figure, fraction=0.9, tightLayout=True):
    """
    Fraction is fraction of height all subplots will be after adjustment.
    Default is 0.9 for 90% of height. 
    """
    # we have to set tight layout before adjustment otherwise if called
    # later it will cause an overlap with the machine diagram
    if (tightLayout):
        _plt.tight_layout()
    
    axs = figure.get_axes()
    
    for ax in axs:
        bbox = ax.get_position()
        bbox.y0 = bbox.y0 * fraction
        bbox.y1 = bbox.y1 * fraction
        ax.set_position(bbox)    

def AddMachineLatticeFromSurveyToFigure(figure, *args, **kwargs):
    """
    kwargs - 'tightLayout' is set to True by default - can be supplied
              in kwargs to force it to false.

    """
    # options
    tightLayout = True
    if 'tightLayout' in kwargs:
        tightLayout = kwargs['tightLayout']

    axoptics  = figure.get_axes()[0]
    _AdjustExistingAxes(figure, tightLayout=tightLayout)
    axmachine = _PrepareMachineAxes(figure)
    
    #concatenate machine lattices
    sf = _CheckItsBDSAsciiData(args[0])
    if len(args) > 1:
        for machine in args[1:]:
            sf.ConcatenateMachine(machine)

    _DrawMachineLattice(axmachine,sf)

    #put callbacks for linked scrolling
    def MachineXlim(ax): 
        axmachine.set_autoscale_on(False)
        axoptics.set_xlim(axmachine.get_xlim())

    def Click(a) : 
        if a.button == 3 : 
            print('Closest element: ',sf.NameFromNearestS(a.xdata))
            
    MachineXlim(axmachine)
    axmachine.callbacks.connect('xlim_changed', MachineXlim)
    figure.canvas.mpl_connect('button_press_event', Click)

def _DrawMachineLattice(axesinstance,bdsasciidataobject):
    ax  = axesinstance #handy shortcut
    bds = bdsasciidataobject

    if not hasattr(bds,"SStart"):
        raise ValueError("This file doesn't have the required column SStart")
    if not hasattr(bds,"ArcLength"):
        raise ValueError("This file doesn't have the required column ArcLength")
    
    def DrawBend(start,length,color='b',alpha=1.0):
        br = _patches.Rectangle((start,-0.1),length,0.2,color=color,alpha=alpha)
        ax.add_patch(br)
    def DrawHKicker(start, length, color='purple', alpha=1.0):
        br = _patches.Rectangle((start,-0.1),length,0.2,color=color,alpha=alpha)
        ax.add_patch(br)
    def DrawVKicker(start, length, color='magenta', alpha=1.0):
        br = _patches.Rectangle((start,-0.1),length,0.2,color=color,alpha=alpha)
        ax.add_patch(br)
    def DrawQuad(start,length,k1l,color='r',alpha=1.0):
        #survey file doesn't have k values
        if k1l > 0 :
            qr = _patches.Rectangle((start,0),length,0.2,color=color,alpha=alpha)
        elif k1l < 0: 
            qr = _patches.Rectangle((start,-0.2),length,0.2,color=color,alpha=alpha)
        else:
            #quadrupole off
            qr = _patches.Rectangle((start,-0.1),length,0.2,color='#B2B2B2',alpha=0.5) #a nice grey in hex
        ax.add_patch(qr)
    def DrawHex(start,length,color,alpha=1.0):
        s = start
        l = length
        edges = _np.array([[s,-0.1],[s,0.1],[s+l/2.,0.13],[s+l,0.1],[s+l,-0.1],[s+l/2.,-0.13]])
        sr = _patches.Polygon(edges,color=color,fill=True,alpha=alpha)
        ax.add_patch(sr)
    def DrawRect(start,length,color,alpha=1.0):
        rect = _patches.Rectangle((start,-0.1),length,0.2,color=color,alpha=alpha)
        ax.add_patch(rect)
    def DrawLine(start,color,alpha=1.0):
        ax.plot([start,start],[-0.2,0.2],'-',color=color,alpha=alpha)
            
    # plot beam line
    smax = bds.SEnd()[-1]
    ax.plot([0,smax],[0,0],'k-',lw=1)
    ax.set_ylim(-0.2,0.2)
 
    # loop over elements and Draw on beamline
    types   = bds.Type()
    lengths = bds.ArcLength()
    starts  = bds.SStart()
    if hasattr(bds,'k1'):
        k1  = bds.k1()
    elif hasattr(bds,'K1'):
        k1  = bds.K1()
    for i in range(len(bds)):
        kw = types[i]
        if kw == 'quadrupole': 
            DrawQuad(starts[i],lengths[i],k1[i], u'#d10000') #red
        elif kw == 'rbend': 
            DrawBend(starts[i],lengths[i], u'#0066cc') #blue
        elif kw == 'sbend': 
            DrawBend(starts[i],lengths[i], u'#0066cc') #blue
        elif kw == 'rcol': 
            DrawRect(starts[i],lengths[i],'k')
        elif kw == 'ecol': 
            DrawRect(starts[i],lengths[i],'k')
        elif kw == 'degrader': 
            DrawRect(starts[i],lengths[i],'k')
        elif kw == 'sextupole':
            DrawHex(starts[i],lengths[i], u'#ffcc00') #yellow
        elif kw == 'octupole':
            DrawHex(starts[i],lengths[i], u'#00994c') #green
        elif kw == 'decapole':
            DrawHex(starts[i],lengths[i], u'#4c33b2') #purple
        elif kw == 'hkick':
            DrawHKicker(starts[i],lengths[i], u'#4c33b2') #purple
        elif kw == 'vkick':
            DrawVKicker(starts[i],lengths[i], u'#ba55d3') #medium orchid
        elif kw == 'drift':
            pass
        elif kw == 'multipole':
            DrawHex(starts[i],lengths[i],'grey',alpha=0.5)
        else:
            #unknown so make light in alpha
            if lengths[i] > 1e-1:
                DrawRect(starts[i],lengths[i],'#cccccc',alpha=0.1) #light grey
            else:
                #relatively short element - just draw a line
                DrawLine(starts[i],'#cccccc',alpha=0.1)


# Predefined lists of tuples for making the standard plots,
# format = (optical_var_name, optical_var_error_name, legend_name)

_BETA = [("Beta_x", "Sigma_Beta_x", r'$\beta_{x}$'),
         ("Beta_y", "Sigma_Beta_y", r'$\beta_{y}$')]

_ALPHA = [("Alpha_x", "Sigma_Alpha_x", r"$\alpha_{x}$"),
          ("Alpha_y", "Sigma_Alpha_y", r"$\alpha_{y}$")]

_DISP = [("Disp_x", "Sigma_Disp_x", r"$D_{x}$"),
         ("Disp_y", "Sigma_Disp_y", r"$D_{y}$")]

_DISP_P = [("Disp_xp", "Sigma_Disp_xp", r"$D_{p_{x}}$"),
           ("Disp_yp", "Sigma_Disp_yp", r"$D_{p_{y}}$")]

_SIGMA = [("Sigma_x", "Sigma_Sigma_x", r"$\sigma_{x}$"),
          ("Sigma_y", "Sigma_Sigma_y", r"$\sigma_{y}$")]

_SIGMA_P = [("Sigma_xp", "Sigma_Sigma_xp", r"$\sigma_{xp}$"),
            ("Sigma_yp", "Sigma_Sigma_yp", r"$\sigma_{yp}$")]

_MEAN = [("Mean_x", "Sigma_Mean_x", r"$\bar{x}$"),
         ("Mean_y", "Sigma_Mean_y", r"$\bar{y}$")]

def _make_plotter(plot_info_tuples, x_label, y_label, title):
    def f_out(bds, outputfilename= None, survey=None, **kwargs):
        sf = _CheckItsBDSAsciiData(bds)

        # options
        tightLayout = True
        if 'tightLayout' in kwargs:
            tightLayout = kwargs['tightLayout']

        # Get the initial N for the two sources
        first_nparticles = sf.Npart()[0]

        plot = _plt.figure(title, **kwargs)
        # Loop over the variables in plot_info_tuples and draw the plots.
        for var, error, legend_name in plot_info_tuples:
            _plt.errorbar(sf.GetColumn('S'),
                          sf.GetColumn(var),
                          yerr=sf.GetColumn(error),
                          label="{}; {}; N = {:.1E}".format(
                              "", legend_name, first_nparticles),
                          capsize=3, **kwargs)

        # Set axis labels and draw legend
        axes = _plt.gcf().gca()
        axes.set_ylabel(y_label)
        axes.set_xlabel(x_label)
        axes.legend(loc='best')

        if survey is not None:
            AddMachineLatticeFromSurveyToFigure(plot, survey)
        if (tightLayout):
            _plt.tight_layout()

        _plt.show(block=False)

        if outputfilename != None:
            if '.' in outputfilename:
                outputfilename = outputfilename.split('.')[0]
            _plt.savefig(outputfilename + '.pdf')
            _plt.savefig(outputfilename + '.png')

        return plot
    return f_out

PlotBeta   = _make_plotter(_BETA,    "S / m", r"$\beta_{x,y}$ / m",      "Beta")
PlotAlpha  = _make_plotter(_ALPHA,   "S / m", r"$\alpha_{x,y}$ / m",     "Alpha")
PlotDisp   = _make_plotter(_DISP,    "S / m", r"$D_{x,y} / m$",          "Dispersion")
PlotDispP  = _make_plotter(_DISP_P,  "S / m", r"$D_{p_{x},p_{y}}$ / m",  "Momentum_Dispersion")
PlotSigma  = _make_plotter(_SIGMA,   "S / m", r"$\sigma_{x,y}$ / m",     "Sigma")
PlotSigmaP = _make_plotter(_SIGMA_P, "S / m", r"$\sigma_{xp,yp}$ / rad", "SigmaP")
PlotMean   = _make_plotter(_MEAN,    "S / m", r"$\bar{x}, \bar{y}$ / m", "Mean")


def PlotBdsimOptics(bdsdata, outputfilename=None, survey=None, **kwargs):
    """
    Display all the optical function plots for a rebdsim optics root file.
    """
    PlotBeta(bdsdata, survey=survey, outputfilename=outputfilename, **kwargs)
    PlotAlpha(bdsdata, survey=survey, outputfilename=outputfilename, **kwargs)
    PlotDisp(bdsdata, survey=survey, outputfilename=outputfilename, **kwargs)
    PlotDispP(bdsdata, survey=survey, outputfilename=outputfilename, **kwargs)
    PlotSigma(bdsdata, survey=survey, outputfilename=outputfilename, **kwargs)
    PlotSigmaP(bdsdata, survey=survey, outputfilename=outputfilename, **kwargs)
    PlotMean(bdsdata, survey=survey, outputfilename=outputfilename, **kwargs)
