import matplotlib.pyplot as plt
import numpy as np


def histogram1D(ax,Data_X,nBins):

    ax.hist(Data_X,nBins, normed=True ,histtype='step')
    ax.grid(True)
    ax.set_xlabel(Data_X.columns)
    ax.set_ylabel('Counts')
