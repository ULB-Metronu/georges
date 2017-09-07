""" Compute some functions for statistics as confidence bounds"""
import numpy as np
from lmfit import Model


def gauss(x, *p):
    a, mu, sigma = p
    return a*np.exp(-(x-mu)**2/(2.*sigma**2))


def gaussian(x, amp, cen, wid):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp/(np.sqrt(2*np.pi)*wid)) * np.exp(-(x-cen)**2 /(2*wid**2))


def gaussian_fit(data, **kwargs):

    if kwargs.get('default_lim',True):
        lim = np.arange(-40, 40, 80/50)

    else:
        lim = kwargs.get('lim')

    y, bin_edges = np.histogram(data, density=True, bins=lim)
    x = (bin_edges[:-1] + bin_edges[1:]) / 2

    gmodel = Model(gaussian)
    result = gmodel.fit(y, x=x, amp=1, cen=x.mean(), wid=x.std())

    if kwargs.get('fit_report', False):
        print(result.fit_report())

    return x, result

#
# def gaussian_fit(data, nbin=50):
#
#     data_bin = np.arange(data.min(), data.max(), (data.max() - data.min()) / nbin)
#     n, bins = np.histogram(data, bins=data_bin, normed=True)
#     bincenter = (bins[:-1] + bins[1:]) / 2
#     coeff, var_matrix = curve_fit(gauss, bincenter, n, p0=[1, data.mean(), data.std()], maxfev=1000)
#     hist_fit = gauss(bincenter, *coeff)
#     error = np.sqrt(np.diag(var_matrix))
#
#     return bincenter, coeff, hist_fit, error

###################################################################

# def mean_confidence_interval(data, confidence=0.95):
#     a = 1.0*np.array(data)
#     n = len(a)
#     m, se = np.mean(a), scipy.stats.sem(a)
#     h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
#     return m, m-h, m+h
#
#
# def std_confidence_interval(data, alpha=0.05):
#
#     n=len(data)
#     s2=np.var(data)
#
#     chi2min=chi2.isf(alpha/2, n-1)
#     chi2max=chi2.isf(1-alpha/2, n-1)
#
#
#     confidence_min=np.sqrt((n-1)*s2/chi2min)
#     confidence_max=np.sqrt((n-1)*s2/chi2max)
#
#     return confidence_min,confidence_max