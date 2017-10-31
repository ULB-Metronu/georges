import matplotlib.pyplot as plt
import numpy as np


def histogram_1dimension(ax,data,nbins,axis_range,datalabel):
    """ Plot an histogram based on panda dataframe """
    if len(axis_range) ==0:
        axis_range=(min(data),max(data))
    ax.hist(data.values, nbins, normed=True, histtype='step', range=axis_range,label=datalabel)
    ax.grid(True)
    ax.set_xlabel(data.name)
    ax.set_ylabel('Counts')