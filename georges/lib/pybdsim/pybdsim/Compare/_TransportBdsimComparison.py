import pymadx as _pymadx
#import pybdsim as _pybdsim
import matplotlib.pyplot as _plt
import numpy as _np

def TransportVsBDSIM(parameter, bdsfile, transfile, transscaling=1, lattice=None, ylabel=None, outputfilename=None):
    #bdsim parameter names and equivalent tfs names
    okParams = ['Sigma_x',
                'Sigma_y',
                'Sigma_xp',
                'Sigma_yp',
                'Beta_x',
                'Beta_y',
                'Disp_x',
                'Disp_y',
                'Alph_x',
                'Alph_y']
    if parameter not in okParams:
        raise ValueError(parameter +' is not a plottable parameter.')

    #check inputs and load data
    bds   = _pybdsim._General.CheckItsBDSAsciiData(bdsim)
    trans = _pymadx.Data.CheckItsTfs(tfs) # TBC is this meant to be MADX Tfs?

    #name of parameter error
    errorparam = 'Sigma_' + _string.lower(parameter)

    maxes = []
    mins  = []
    #set plot minimum to 0 for any parameter which can't be negative.
    if parameter[:4] != ('Disp' or 'Alph'):
            mins.append(0)

    _plt.figure()
    if parameter in bds.names:
        data   = bds.GetColumn(parameter)
        bdsimlabel = 'BDSIM, NPrimaries = %1.0e' %bds.GetColumn('Npart')[0]
        #Check if errors exist
        if errorparam in bds.names:
            errors = bds.GetColumn(errorparam)
            _plt.errorbar(bds.GetColumn('S'), data, yerr=errors, label=bdsimlabel)
        else:
            _plt.plot(bds.GetColumn('S'), data,label=bdsimlabel)
        maxes.append(_np.max(data + errors))
        mins.append(_np.max(data - errors))
    if parameter in trans.names:
        data = trans.GetColumn(parameter)*transscaling
        _plt.plot(trans.GetColumn('S'), data, label='TRANSPORT')
        maxes.append(_np.max(data))
        mins.append(_np.min(data))

    _plt.xlabel('S (m)')
    if ylabel != None:
        _plt.ylabel(ylabel)
    _plt.legend(loc=0)
    _plt.ylim(1.1*_np.min(mins), 1.1*_np.max(maxes))
    if lattice != None:
        _pybdsim.Plot.AddMachineLatticeFromSurveyToFigure(_plt.gcf(),lattice)
    if outputfilename != None:
        if '.' in outputfilename:
            outputfilename = outputfilename.split('.')[0]
        _plt.savefig(outputfilename+'.pdf')
        _plt.savefig(outputfilename+'.png')
