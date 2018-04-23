"""
Functions for comparing the optical distributions of two
BDSIM models.

Functions for plotting the individual optical functions, and an
eighth, helper function ''compare_all_optics``, to plot display all
seven in one go.
"""

import os.path as _path
import matplotlib.pyplot as _plt

from .. import Data
from .. import Plot

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

# use closure to avoid tonnes of boilerplate code as happened with
# MadxBdsimComparison.py
def _make_plotter(plot_info_tuples, x_label, y_label, title):
    def f_out(first, second, first_name=None, second_name=None,
              survey=None, **kwargs):
        """first and second should be rebdsimOptics output root files."""
        if not _path.isfile(first):
            raise IOError("file \"{}\" not found!".format(first))
        if not _path.isfile(second):
            raise IOError("file \"{}\" not found!".format(second))

        # If no names provided then just use the filenames.
        first_name = (_path.splitext(_path.basename(first))[0]
                      if first_name is None else first_name)
        second_name = (_path.splitext(_path.basename(second))[0]
                       if second_name is None else second_name)

        first = pybdsim.Data.Load(first)
        second = pybdsim.Data.Load(second)

        # Get the initial N for the two sources
        first_nparticles = first.Npart()[0]
        second_nparticles = second.Npart()[0]

        plot = _plt.figure(title, **kwargs)
        # Loop over the variables in plot_info_tuples and draw the plots.
        for var, error, legend_name in plot_info_tuples:
            _plt.errorbar(first.GetColumn('S'),
                          first.GetColumn(var),
                          yerr=first.GetColumn(error),
                          label="{}; {}; N = {:.1E}".format(
                              first_name, legend_name, first_nparticles),
                          capsize=3, **kwargs)
            _plt.errorbar(second.GetColumn('S'),
                          second.GetColumn(var),
                          yerr=second.GetColumn(error),
                          label="{}; {}; N = {:.1E}".format(
                              second_name, legend_name, second_nparticles),
                          capsize=3, **kwargs)

        # Set axis labels and draw legend
        axes = _plt.gcf().gca()
        axes.set_ylabel(y_label)
        axes.set_xlabel(x_label)
        axes.legend(loc='best')

        if survey is not None:
            pybdsim.Plot.AddMachineLatticeFromSurveyToFigure(plot, survey)
        _plt.show(block=False)
        return plot
    return f_out

PlotBeta   = _make_plotter(_BETA,    "S / m", r"$\beta_{x,y}$ / m",      "Beta")
PlotAlpha  = _make_plotter(_ALPHA,   "S / m", r"$\alpha_{x,y}$ / m",     "Alpha")
PlotDisp   = _make_plotter(_DISP,    "S / m", r"$D_{x,y} / m$",          "Dispersion")
PlotDispP  = _make_plotter(_DISP_P,  "S / m", r"$D_{p_{x},p_{y}}$ / m",  "Momentum_Dispersion")
PlotSigma  = _make_plotter(_SIGMA,   "S / m", r"$\sigma_{x,y}$ / m",     "Sigma")
PlotSigmaP = _make_plotter(_SIGMA_P, "S / m", r"$\sigma_{xp,yp}$ / rad", "SigmaP")
PlotMean   = _make_plotter(_MEAN,    "S / m", r"$\bar{x}, \bar{y}$ / m", "Mean")


def BDSIMVsBDSIM(first, second, first_name=None,
                 second_name=None, survey=None, **kwargs):
    """
    Display all the optical function plots for the two input optics files.
    """
    PlotBeta(first, second, first_name=first_name,
             second_name=second_name, survey=survey, **kwargs)
    PlotAlpha(first, second, first_name=first_name,
              second_name=second_name, survey=survey, **kwargs)
    PlotDisp(first, second, first_name=first_name,
             second_name=second_name, survey=survey, **kwargs)
    PlotDispP(first, second, first_name=first_name,
              second_name=second_name, survey=survey, **kwargs)
    PlotSigma(first, second, first_name=first_name,
              second_name=second_name, survey=survey, **kwargs)
    PlotSigmaP(first, second, first_name=first_name,
               second_name=second_name, survey=survey, **kwargs)
    PlotMean(first, second, first_name=first_name,
             second_name=second_name, survey=survey, **kwargs)
