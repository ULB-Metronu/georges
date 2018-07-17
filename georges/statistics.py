""" Compute some functions for statistics as confidence bounds"""
import numpy as np
from lmfit.models import GaussianModel
from lmfit.models import LorentzianModel


REGISTERED_MODEL = ['Gaussian', 'Lorentzian']


class StatisticException(Exception):
    """Exception raised for errors in the beam plotting module."""

    def __init__(self, m):
        self.message = m


def histogram_fit(data, bounds_binning=50, verbose=False, **kwargs):

    """ All models are available on https://lmfit.github.io/lmfit-py/builtin_models.html#lmfit.models"""
    y, bin_edges = np.histogram(data, density=False, bins=bounds_binning)
    x = (bin_edges[:-1] + bin_edges[1:]) / 2

    model = kwargs.get("model", "Gaussian")

    if model not in REGISTERED_MODEL:
        # TODO make an excception
        raise StatisticException(f"The model {model} is not in {REGISTERED_MODEL}")

    if model == "Gaussian":
        mod = GaussianModel()

    if model == "Lorentzian":
        mod = LorentzianModel()

    result = mod.fit(data=y, x=x, center=data.mean(), sigma=data.var())

    if verbose:
        print(result.fit_report())
    return x, result
