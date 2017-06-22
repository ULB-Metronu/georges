""" Compute some functions for statistics as confidence bounds"""
import numpy as np
import scipy as sp
import scipy.stats
from scipy.stats import chi2

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

def std_confidence_interval(data, alpha=0.05):
    
    n=len(data)
    s2=np.var(data)
    
    chi2min=chi2.isf(alpha/2, n-1)
    chi2max=chi2.isf(1-alpha/2, n-1)
    
    
    confidence_min=np.sqrt((n-1)*s2/chi2min)
    confidence_max=np.sqrt((n-1)*s2/chi2max)
    
    return confidence_min,confidence_max