import pymadx as _pymadx
#from .. import pybdsim as _pybdsim
import matplotlib.pyplot as _plt
import numpy as _np
from os.path import isfile
from matplotlib.backends.backend_pdf import PdfPages
import datetime

def MadxVsBDSIM(tfs, bdsim, survey=None, functions=None,
                postfunctions=None, figsize=(12, 5), saveAll=True):
    """
    Compares MadX and BDSIM optics variables.
    User must provide a tfsoptIn file or Tfsinstance and a BDSAscii file or instance.

    +-----------------+---------------------------------------------------------+
    | **Parameters**  | **Description**                                         |
    +-----------------+---------------------------------------------------------+
    | tfs             | Tfs file or pymadx.Data.Tfs instance.                   |
    +-----------------+---------------------------------------------------------+
    | bdsim           | Optics root file (from rebdsimOptics or rebdsim).       |
    +-----------------+---------------------------------------------------------+
    | survey          | BDSIM model survey.                                     |
    +-----------------+---------------------------------------------------------+
    | functions       | Hook for users to add their functions that are called   |
    |                 | immediately prior to the addition of the plot. Use a    |
    |                 | lambda function to add functions with arguments. Can    |
    |                 | be a function or a list of functions.                   |
    +-----------------+---------------------------------------------------------+
    | figsize         | Figure size for all figures - default is (12,5)         |
    +-----------------+---------------------------------------------------------+
    
    """

    _CheckFilesExist(tfs, bdsim, survey)

    tfsinst   = _pymadx.Data.CheckItsTfs(tfs)
    bdsinst   = _pybdsim._General.CheckItsBDSAsciiData(bdsim)

    tfsopt    = _GetTfsOptics(tfsinst)
    bdsopt    = _GetBDSIMOptics(bdsinst)

    if survey is None:
        survey = tfsinst

    if saveAll:
        # A4 size with quarter inch margin on sides and on top to
        # ensure figures don't come out looking cramped and bad.
        figsize=(11.44, 8.02)

    figures = [PlotBetas(tfsopt, bdsopt, survey=survey,
                         functions=functions,
                         postfunctions=postfunctions,
                         figsize=figsize),
               PlotAlphas(tfsopt, bdsopt, survey=survey,
                          functions=functions,
                          postfunctions=postfunctions,
                          figsize=figsize),
               PlotDs(tfsopt, bdsopt, survey=survey,
                      functions=functions,
                      postfunctions=postfunctions,
                      figsize=figsize),
               PlotDps(tfsopt, bdsopt, survey=survey,
                       functions=functions,
                       postfunctions=postfunctions,
                       figsize=figsize),
               PlotSigmas(tfsopt, bdsopt, survey=survey,
                          functions=functions,
                          postfunctions=postfunctions,
                          figsize=figsize),
               PlotSigmasP(tfsopt, bdsopt, survey=survey,
                           functions=functions,
                           postfunctions=postfunctions,
                           figsize=figsize),
               PlotMeans(tfsopt, bdsopt, survey=survey,
                         functions=functions,
                         postfunctions=postfunctions,
                         figsize=figsize)]

    if saveAll:
        tfsname = repr(tfsinst)
        bdsname = repr(bdsinst)
        output_filename = "optics-report.pdf"
        # Should have a more descriptive name really.
        with PdfPages(output_filename) as pdf:
            for figure in figures:
                pdf.savefig(figure)
            d = pdf.infodict()
            d['Title'] = "{} (TFS) VS {} (BDSIM) Optical Comparison".format(
                tfsname, bdsname)
            d['CreationDate'] = datetime.datetime.today()

        print("Written ", output_filename)

def MadxVsBDSIMOrbit(tfs, bdsim, survey=None, functions=None, postfunctions=None):

    _CheckFilesExist(tfs, bdsim, survey)
    tfsinst   = _pymadx.Data.CheckItsTfs(tfs)
    bdsinst   = _pybdsim._General.CheckItsBDSAsciiData(bdsim)

    tfsorbit  = _GetTfsOrbit(tfsinst)
    bdsopt    = _GetBDSIMOptics(bdsinst)

    if survey is None:
        survey = tfsinst

    PlotOrbit(tfsorbit, bdsopt, survey=survey, functions=functions)
    #PlotOrbitResiduals(tfsorbit, bdsopt, survey=survey, functions=functions)

def PrepareResiduals(tfs, bds, survey=None, verbose=False):
    """
    Filter the tfs to provide data that will match the 
    BDSIM data in s position. 
    """
    _CheckFilesExist(tfs, bds, survey)
    tfsinst   = _pymadx.Data.CheckItsTfs(tfs)
    bdsinst   = _pybdsim._General.CheckItsBDSAsciiData(bds) # works for root files too

    bdss = bdsinst.s()

    keys = ['S', 'X', 'PX', 'Y', 'PY']

    tfsdata = {
        'S':[],
        'X':[],
        'PX':[],
        'Y':[],
        'PY':[]
        }
    
    for s in bdss:
        ind = tfsinst.IndexFromNearestS(s)
        if (verbose):
            print('bdsim s:',s,' index in tfs:',ind, tfsinst[ind])
        for key in keys:
            tfsdata[key].append(tfsinst[ind][key])

    for key in keys:
        tfsdata[key] = _np.array(tfsdata[key])

    return tfsdata

    
def MadxVsBDSIMFromGMAD(tfs, gmad):
    """Runs the BDSIM model provided by the gmad file given, gets the
    optics, and then compares them with TFS.

    tfs - path to TFS file or a pymadx.Data.Tfs instance
    gmad - path to gmad file to be run.

    """

    bdsimOptics = _pybdsim.Run.GetOpticsFromGMAD(gmad)
    # Use TFS for survey but perhaps one day can directly get the
    # model from the ROOT output (having modified the above
    # function..  Currently there's no interface for this.
    MadxVsBDSIM(tfs, bdsimOptics, survey=tfs)


def _GetBDSIMOptics(optics):
    '''
    Takes a BDSAscii instance.
    Return a dictionary of lists matching the variable with the list of values.
    '''
    optvars = {}
    for variable in optics.names:
        datum = getattr(optics, variable)()
        optvars[variable] = datum
    return optvars

def _GetTfsOptics(optics):
    '''
    Takes either Tfs instance.  Returns dictionary of lists.
    '''

    MADXOpticsVariables = frozenset(['NAME',
                                     'S',
                                     'BETX',
                                     'BETY',
                                     'ALFX',
                                     'ALFY',
                                     'DX',
                                     'DPX',
                                     'DY',
                                     'DPY',
                                     'SIGMAX',
                                     'SIGMAY',
                                     'SIGMAXP',
                                     'SIGMAYP',
                                     'X',
                                     'Y',
                                     'PX',
                                     'PY'])

    optvars = {}
    for variable in MADXOpticsVariables:
        optvars[variable] = optics.GetColumn(variable)
    return optvars

def _GetTfsOrbit(optics):
    '''
    Takes either Tfs instance.  Returns dictionary of lists.
    '''

    MADXOpticsVariables = frozenset(['S',
                                     'X',
                                     'Y',
                                     'PX',
                                     'PY'])

    optvars = {}
    for variable in MADXOpticsVariables:
        optvars[variable] = optics.GetColumn(variable)
    return optvars

def PlotBetas(tfsopt, bdsopt, survey=None, functions=None, postfunctions=None, figsize=(12,5)):

    N = str(int(bdsopt['Npart'][0]))  #number of primaries.
    betaPlot = _plt.figure('Beta', figsize=figsize)
    _plt.plot(tfsopt['S'], tfsopt['BETX'], 'b', label=r'MADX $\beta_{x}$')
    _plt.plot(tfsopt['S'], tfsopt['BETY'], 'g', label=r'MADX $\beta_{y}$')

    #bds
    _plt.errorbar(bdsopt['S'], bdsopt['Beta_x'],
                  yerr=bdsopt['Sigma_Beta_x'],
                  label=r'BDSIM $\beta_{x}$' + ' ; N = ' + N,
                  marker='x',
                  ls = '',
                  color='b')

    _plt.errorbar(bdsopt['S'], bdsopt['Beta_y'],
                  yerr=bdsopt['Sigma_Beta_y'],
                  label=r'BDSIM $\beta_{y}$' + ' ; N = ' + N,
                  marker='x',
                  ls = '',
                  color='g')

    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$\beta_{x,y}$ / m')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')

    print("survey =", survey)

    _CallUserFigureFunctions(functions)
    _AddSurvey(betaPlot, survey)
    _CallUserFigureFunctions(postfunctions)
    
    _plt.show(block=False)
    return betaPlot

def PlotAlphas(tfsopt, bdsopt, survey=None, functions=None, postfunctions=None, figsize=(12,5)):
    N = str(int(bdsopt['Npart'][0]))  #number of primaries.
    alphaPlot = _plt.figure('Alpha', figsize=figsize)
    #tfs
    _plt.plot(tfsopt['S'], tfsopt['ALFX'], 'b', label=r'MADX $\alpha_{x}$')
    _plt.plot(tfsopt['S'], tfsopt['ALFY'], 'g', label=r'MADX $\alpha_{y}$')
    #bds
    _plt.errorbar(bdsopt['S'], bdsopt['Alpha_x'],
                  yerr=bdsopt['Sigma_Alpha_x'],
                  label=r'BDSIM $\alpha_{x}$' + ' ; N = ' + N,
                  fmt='b.', capsize=3)

    _plt.errorbar(bdsopt['S'], bdsopt['Alpha_y'],
                  yerr=bdsopt['Sigma_Alpha_y'],
                  label=r'BDSIM $\alpha_{y}$' + ' ; N = ' + N,
                  fmt='g.', capsize=3)

    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$\alpha_{x,y}$ / m')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')

    _CallUserFigureFunctions(functions)
    _AddSurvey(alphaPlot, survey)
    _CallUserFigureFunctions(postfunctions)
    
    _plt.show(block=False)
    return alphaPlot

def PlotDs(tfsopt, bdsopt, survey=None, functions=None, postfunctions=None, figsize=(12,5)):
    N = str(int(bdsopt['Npart'][0]))  #number of primaries.
    dispPlot = _plt.figure('Dispersion', figsize=figsize)
    #tfs
    _plt.plot(tfsopt['S'], tfsopt['DX'], 'b', label=r'MADX $D_{x}$')
    _plt.plot(tfsopt['S'], tfsopt['DY'], 'g', label=r'MADX $D_{y}$')
    #bds
    _plt.errorbar(bdsopt['S'], bdsopt['Disp_x'],
                  yerr=bdsopt['Sigma_Disp_x'],
                  label=r'BDSIM $D_{x}$' + ' ; N = ' + N,
                  fmt='b.', capsize=3)

    _plt.errorbar(bdsopt['S'], bdsopt['Disp_y'],
                  yerr=bdsopt['Sigma_Disp_y'],
                  label=r'BDSIM $D_{y}$' + ' ; N = ' + N,
                  fmt='g.', capsize=3)

    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$D_{x,y} / m$')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')

    _CallUserFigureFunctions(functions)
    _AddSurvey(dispPlot, survey)
    _CallUserFigureFunctions(postfunctions)
    
    _plt.show(block=False)
    return dispPlot

def PlotDps(tfsopt, bdsopt, survey=None, functions=None, postfunctions=None, figsize=(12,5)):
    N = str(int(bdsopt['Npart'][0]))  #number of primaries.
    dispPPlot = _plt.figure('Momentum_Dispersion', figsize=figsize)
    #tfs
    _plt.plot(tfsopt['S'], tfsopt['DPX'], 'b', label=r'MADX $D_{p_{x}}$')
    _plt.plot(tfsopt['S'], tfsopt['DPY'], 'g', label=r'MADX $D_{p_{y}}$')
    #bds
    _plt.errorbar(bdsopt['S'], bdsopt['Disp_xp'],
                  yerr=bdsopt['Sigma_Disp_xp'],
                  label=r'BDSIM $D_{p_{x}}$' + ' ; N = ' + N,
                  fmt='b.', capsize=3)

    _plt.errorbar(bdsopt['S'], bdsopt['Disp_yp'],
                  yerr=bdsopt['Sigma_Disp_yp'],
                  label=r'BDSIM $D_{p_{y}}$' + ' ; N = ' + N,
                  fmt='g.', capsize=3)

    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$D_{p_{x},p_{y}}$ / m')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')

    _CallUserFigureFunctions(functions)
    _AddSurvey(dispPPlot, survey)
    _CallUserFigureFunctions(postfunctions)
    
    _plt.show(block=False)
    return dispPPlot

def PlotSigmas(tfsopt, bdsopt, survey=None, functions=None, postfunctions=None, figsize=(12,5)):
    N = str(int(bdsopt['Npart'][0]))  #number of primaries.
    sigmaPlot = _plt.figure('Sigma', figsize=figsize)
    #tfs
    _plt.plot(tfsopt['S'], tfsopt['SIGMAX'], 'b', label=r'MADX $\sigma_{x}$')
    _plt.plot(tfsopt['S'], tfsopt['SIGMAY'], 'g', label=r'MADX $\sigma_{y}$')
    #bds
    _plt.errorbar(bdsopt['S'],
                  bdsopt['Sigma_x'],
                  yerr=bdsopt['Sigma_Sigma_x'],
                  label=r'BDSIM $\sigma_{x}$' + ' ; N = ' + N,
                  fmt='b.', capsize=3)

    _plt.errorbar(bdsopt['S'],
                  bdsopt['Sigma_y'],
                  yerr=bdsopt['Sigma_Sigma_y'],
                  label=r'BDSIM $\sigma_{y}$' + ' ; N = ' + N,
                  fmt='g.', capsize=3)
    
    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$\sigma_{x,y}$ / m')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')

    _CallUserFigureFunctions(functions)
    _AddSurvey(sigmaPlot, survey)
    _CallUserFigureFunctions(postfunctions)
    
    _plt.show(block=False)
    return sigmaPlot

def PlotSigmasP(tfsopt, bdsopt, survey=None, functions=None, postfunctions=None, figsize=(12,5)):
    N = str(int(bdsopt['Npart'][0]))  #number of primaries.
    sigmaPPlot = _plt.figure('SigmaP', figsize=figsize)
    #tfs
    _plt.plot(tfsopt['S'], tfsopt['SIGMAXP'], 'b', label=r'MADX $\sigma_{xp}$')
    _plt.plot(tfsopt['S'], tfsopt['SIGMAYP'], 'g', label=r'MADX $\sigma_{yp}$')
    #bds
    _plt.errorbar(bdsopt['S'],
                  bdsopt['Sigma_xp'],
                  yerr=bdsopt['Sigma_Sigma_xp'],
                  label=r'BDSIM $\sigma_{xp}$' + ' ; N = ' + N,
                  fmt='b.', capsize=3)

    _plt.errorbar(bdsopt['S'],
                  bdsopt['Sigma_yp'],
                  yerr=bdsopt['Sigma_Sigma_yp'],
                  label=r'BDSIM $\sigma_{yp}$' + ' ; N = ' + N,
                  fmt='g.', capsize=3)
    
    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$\sigma_{xp,yp}$ / rad')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')

    _CallUserFigureFunctions(functions)
    _AddSurvey(sigmaPPlot, survey)
    _CallUserFigureFunctions(postfunctions)

    _plt.show(block=False)
    return sigmaPPlot

def PlotMeans(tfsopt, bdsopt, survey=None, functions=None, postfunctions=None, figsize=(12,5)):
    N = str(int(bdsopt['Npart'][0]))  #number of primaries.
    meanPlot = _plt.figure('Mean', figsize=figsize)
    #tfs
    _plt.plot(tfsopt['S'], tfsopt['X'], 'b', label=r'MADX $\bar{x}$')
    _plt.plot(tfsopt['S'], tfsopt['Y'], 'g', label=r'MADX $\bar{y}$')

    #bdsim
    _plt.errorbar(bdsopt['S'], bdsopt['Mean_x'],
                  yerr=bdsopt['Sigma_Mean_x'],
                  label=r'BDSIM $\bar{x}$' + ' ; N = ' + N,
                  fmt='b.', capsize=3)

    _plt.errorbar(bdsopt['S'], bdsopt['Mean_y'],
                  yerr=bdsopt['Sigma_Mean_y'],
                  label=r'BDSIM $\bar{y}$' + ' ; N = ' + N,
                  fmt='g.', capsize=3)

    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$\bar{x}, \bar{y}$ / m')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')

    _CallUserFigureFunctions(functions)
    _AddSurvey(meanPlot, survey)
    _CallUserFigureFunctions(postfunctions)
    
    _plt.show(block=False)
    return meanPlot

def PlotOrbit(tfsopt, bdsopt, survey=None, functions=None, postfunctions=None, figsize=(12,5)):
    orbitPlot = _plt.figure('Orbit', figsize=figsize)

    #tfs
    _plt.plot(tfsopt['S'], tfsopt['X'], 'b', label=r'MADX $\bar{x}$')
    _plt.plot(tfsopt['S'], tfsopt['Y'], 'g', label=r'MADX $\bar{y}$')

    #bdsim
    _plt.plot(bdsopt['s'], bdsopt['x'], 'b.', label='BDSIM x')
    _plt.plot(bdsopt['s'], bdsopt['x'], 'b-', alpha=0.4)
    _plt.plot(bdsopt['s'], bdsopt['y'], 'g.', label='BDSIM y')
    _plt.plot(bdsopt['s'], bdsopt['y'], 'g-', alpha=0.4)
    
    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$\bar{x}, \bar{y}$ / m')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')

    _CallUserFigureFunctions(functions)
    _AddSurvey(orbitPlot, survey)
    _CallUserFigureFunctions(postfunctions)

    _plt.show(block=False)
    return orbitPlot

def PlotOrbitResiduals(tfs, bds, survey=None, functions=None, postfunctions=None, verbose=False, figsize=(12,5)):
    _CheckFilesExist(tfs, bds, survey)
    tfsinst   = _pymadx.Data.CheckItsTfs(tfs)
    bdsinst   = _pybdsim._General.CheckItsBDSAsciiData(bds) # works for root files too
    tfsd = PrepareResiduals(tfs, bds)

    if survey is None:
        survey = tfsinst
    
    s   = bdsinst.s()
    dx  = tfsd['X']  - bdsinst.x()
    dxp = tfsd['PX'] - bdsinst.xp()
    dy  = tfsd['Y']  - bdsinst.y()
    dyp = tfsd['PY'] - bdsinst.yp()

    orbRes = _plt.figure('OrbitResiduals', figsize=figsize)
    _plt.plot(s, abs(dx),  '.', label='x')
    _plt.plot(s, abs(dxp), '.', label='xp')
    _plt.plot(s, abs(dy),  '.', label='y')
    _plt.plot(s, abs(dyp), '.', label='yp')
    _plt.legend()
    _plt.xlabel('S (m)')
    _plt.ylabel('|residual| (m, rad)')
    _plt.yscale('log')

    _CallUserFigureFunctions(functions)
    _AddSurvey(orbRes, survey)

    _plt.show(block=False)
    return orbRes

def _AddSurvey(figure, survey):
    if survey is None:
        return
    if isinstance(survey, basestring): # If BDSIM ASCII survey file
        if survey.split(".")[-1] == 'dat':
            _pybdsim.Plot.AddMachineLatticeFromSurveyToFigure(figure,survey)
    # If BDSIM ASCII survey instance
    elif isinstance(survey, _pybdsim.Data.BDSAsciiData):
        _pybdsim.Plot.AddMachineLatticeToFigure(figure,survey)
    elif isinstance(survey, _pymadx.Data.Tfs): # If TFS
        _pymadx.Plot.AddMachineLatticeToFigure(figure,survey)
    # if a (BDSIM) ROOT file
    elif pybdsim._General.IsROOTFile(survey):
        pass

def _ProcessInput(tfsOptics, bdsimOptics):

    if not isinstance(tfsOptics, (_pymadx.Data.Tfs, basestring)):
        raise TypeError("tfsOptics should be either a path to a tfs file or "
                        "a pymadx.Data.Tfs instance!")
    if not isinstance(bdsimOptics, _pybdsim.Data.BDSAsciiData):
        raise TypeError("bdsimOptics should be either be a path to a "
                        "BDSAsciiData file or a pybdsim.Data.BDSAsciiData "
                        "instance")

    if isinstance(tfsOptics, basestring):
        tfsOptics = _pymadx.Data.Tfs(tfsOptics)
    if isinstance(tfsOptics, basestring):
        bdsimOptics = _pybdsim.Data.Load(bdsimOptics)

    return tfsOptics, bdsimOptics

def _CheckFilesExist(tfs, bdsim, survey):
    '''
    Otherwise such errors are too cryptic.
    '''
    if isinstance(tfs, basestring):
        if not isfile(tfs):
            raise IOError("File not found: ", tfs)
    if isinstance(bdsim, basestring) and not isfile(bdsim):
        raise IOError("File not found: ", bdsim)
    if isinstance(survey, basestring) and not isfile(survey):
        raise IOError("File not found: ", survey)


def _CallUserFigureFunctions(functions):
    if isinstance(functions, list):
        for function in functions:
            if callable(function):
                function()
    elif callable(functions):
        functions()
