import numpy as np
from lmfit.models import GaussianModel


class StatisticException(Exception):
    """Exception raised for errors in the beam plotting module."""
    def __init__(self, m):
        self.message = m


def histogram_fit(data, bounds_binning=50, verbose=False, model=GaussianModel):
    """ All models are available on https://lmfit.github.io/lmfit-py/builtin_models.html#lmfit.models"""
    y, bin_edges = np.histogram(data, density=False, bins=bounds_binning)
    x = (bin_edges[:-1] + bin_edges[1:]) / 2
    result = model().fit(data=y, x=x, center=data.mean(), sigma=data.std())
    if verbose:
        print(result.fit_report())
    return x, result
