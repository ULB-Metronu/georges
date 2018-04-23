import pysad as _pysad
import pybdsim as _pybdsim
import matplotlib.pyplot as _plt

def SadComparison(sad, bdsim, survey=None, functions=None, figsize=(12,5)):
    sad = _pysad.Reader.Flat(sad)
    bds = _pybdsim.Data.Load(bdsim)

    PlotAlphas(sad,bds)
    PlotDs(sad,bds)
    PlotDps(sad,bds)
    PlotSigmas(sad,bds)

def PlotSigmas(sad, bds, figsize=(12,5)):
    betaPlot = _plt.figure('Sigma', figsize=figsize)

    # sad
    _plt.plot(sad.dataDict['POSITION'], sad.dataDict['SIGMAX'], 'b', label=r'SAD $\sigma_{x}$')
    _plt.plot(sad.dataDict['POSITION'], sad.dataDict['SIGMAY'], 'g', label=r'SAD $\sigma_{y}$')

    # bds
    _plt.errorbar(bds.S(),
                  bds.Sigma_x(),
                  yerr=bds.Sigma_Sigma_x(),
                  label=r'BDSIM $\sigma_{x}$',
                  fmt='b.', capsize=3)

    # bds
    _plt.errorbar(bds.S(),
                  bds.Sigma_y(),
                  yerr=bds.Sigma_Sigma_y(),
                  label=r'BDSIM $\sigma_{y}$',
                  fmt='g.', capsize=3)


    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$\sigma_{x,y}$ / m')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')

def PlotDs(sad, bds, figsize=(12, 5)):
    dispPlot = _plt.figure('Dispersion', figsize=figsize)

    # sad
    _plt.plot(sad.dataDict['POSITION'], sad.dataDict['EX'], 'b', label=r'SAD $\eta_{x}$')
    _plt.plot(sad.dataDict['POSITION'], sad.dataDict['EY'], 'g', label=r'SAD $\eta_{y}$')

    # bds
    _plt.errorbar(bds.S(),
                  bds.Disp_x(),
                  yerr=bds.Sigma_Disp_x(),
                  label=r'BDSIM $\eta_{x}$',
                  fmt='b.', capsize=3)

    # bds
    _plt.errorbar(bds.S(),
                  bds.Sigma_y(),
                  yerr=bds.Sigma_Disp_y(),
                  label=r'BDSIM $\eta_{y}$',
                  fmt='g.', capsize=3)


    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$\eta_{x,y}$ / m')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')

def PlotDps(sad, bds, figsize=(12, 5)):
    dpsPlot = _plt.figure('Angular dispersion', figsize=figsize)

    # sad
    _plt.plot(sad.dataDict['POSITION'], sad.dataDict['EPX'], 'b', label=r'SAD $\eta_{x}$')
    _plt.plot(sad.dataDict['POSITION'], sad.dataDict['EPY'], 'g', label=r'SAD $\eta_{y}$')

    # bds
    _plt.errorbar(bds.S(),
                  bds.Disp_xp(),
                  yerr=bds.Sigma_Disp_xp(),
                  label=r'BDSIM $\eta_{x}^{\prime}$',
                  fmt='b.', capsize=3)

    # bds
    _plt.errorbar(bds.S(),
                  bds.Sigma_yp(),
                  yerr=bds.Sigma_Disp_yp(),
                  label=r'BDSIM $\eta_{y}^{\prime}$',
                  fmt='g.', capsize=3)


    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$\eta_{x,y}^{\prime}$ / m')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')


def PlotAlphas(sad, bds, figsize=(12, 5)):
    dispPlot = _plt.figure('Alpha', figsize=figsize)

    # sad
    _plt.plot(sad.dataDict['POSITION'], sad.dataDict['AX'], 'b', label=r'SAD $\alpha_{x}$')
    _plt.plot(sad.dataDict['POSITION'], sad.dataDict['AY'], 'g', label=r'SAD $\alpha_{y}$')

    # bds
    _plt.errorbar(bds.S(),
                  bds.Alpha_x(),
                  yerr=bds.Sigma_Alpha_x(),
                  label=r'BDSIM $\alpha_{x}$',
                  fmt='b.', capsize=3)

    # bds
    _plt.errorbar(bds.S(),
                  bds.Alpha_y(),
                  yerr=bds.Sigma_Alpha_y(),
                  label=r'BDSIM $\alpha_{y}$',
                  fmt='g.', capsize=3)


    axes = _plt.gcf().gca()
    axes.set_ylabel(r'$\alpha_{x,y}$ / m')
    axes.set_xlabel('S / m')
    axes.legend(loc='best')
