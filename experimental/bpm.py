import re, os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


DATA_BPM_PCVUE = {
    'max_current': {
        'regex': re.compile("Max. current : (.*) nA/mm"),
    },
    'average_integrated_current': {
        'regex': re.compile("Average Integrated Current : (.*) nA"),
    },
    'centroid': {
        'regex': re.compile("centroid \(mm\) :(.*)"),
    },
    'sigma': {
        'regex': re.compile("sigma \(mm\) :(.*)"),
    },
    'chisqr': {
        'regex': re.compile("chisqr :(.*)"),
    },
    'channel': {
        'regex': re.compile("(..)		TRUE		(.*)		(.*)"),
    },
}


def read_data_file(file, regex_data, flag, path='', dtype=float):
    data = {}
    for k in regex_data:
        data[k] = {'data': {}, 'regex': regex_data[k]['regex']}
        v = data[k]
        if flag.get('init'):
            v['data'][flag['init']] = []
        if flag.get('next'):
            v['data'][flag['next']] = []
        f = flag.get('init', 'init')
        for line in open(os.path.join(path, file)):
            if line.startswith(flag['trigger']):
                f = flag.get('next', 'next')
            for match in re.finditer(v['regex'], line):
                d = list(map(str.strip, match.groups()))
                v['data'][f].append(list(map(dtype,d)))
        v['data'][flag.get('init', 'init')] = np.array(v['data'].get(flag.get('init', 'init')))
        v['data'][flag.get('next', 'next')] = np.array(v['data'].get(flag.get('next', 'next')))
    return data


def gaussian(x, a, mu, sigma):
    return a * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


def bpm_fit(x, y):
    ar = np.trapz(y / np.sum(y) * len(y), x)
    mean = np.mean(x * y / np.sum(y) * len(y))
    rms = np.std(x * y / np.sum(y) * len(y))

    popt, pcov = curve_fit(gaussian, x, y, p0=[ar, mean, rms])

    return [popt, np.sqrt(pcov.diagonal())]


def bpm(**kwargs):
    """Process a PCVue BPM data file and read associated data onto a DataFrame.
    :param kwargs: parameters are:
        - path: the path were the data file and bpm data files can be found
        - file: the data file name
        - bpm: the BPM name
        - instrument: the specific parsing data for the instrument
    """
    path = kwargs.get("path", ".")
    file = kwargs.get("file")
    bpm_name = kwargs.get("bpm", "BPM")
    instrument = kwargs.get("instrument", DATA_BPM_PCVUE)
    try:
        bpm_data = pd.read_csv(os.path.join(path, file), skiprows=0, encoding="utf-8-sig")
    except:
        bpm_data = pd.read_excel(os.path.join(path, file), skiprows=0, encoding="utf-8-sig")
    bpm_data["{}_data".format(bpm_name)] = bpm_data[bpm_name].apply(
        lambda x: read_data_file(
            "{0:04d}.txt".format(x),
            instrument,
            {'init': 'Y', 'next': 'X', 'trigger': "Y-Wire data"},
            path
        )
    )
    bpm_data["{}_fit_X".format(bpm_name)] = bpm_data["{}_data".format(bpm_name)].apply(
        lambda x: bpm_fit(
            x['channel']['data']['X'][:, 1],
            x['channel']['data']['X'][:, 2]
        )
    )
    bpm_data["{}_fit_Y".format(bpm_name)] = bpm_data["{}_data".format(bpm_name)].apply(
        lambda x: bpm_fit(
            x['channel']['data']['Y'][:, 1],
            x['channel']['data']['Y'][:, 2]
        )
    )
    bpm_data["{}_fit_sigma_X".format(bpm_name)] = bpm_data["{}_fit_X".format(bpm_name)].apply(lambda x: np.abs(x[0][2]))
    bpm_data["{}_fit_sigma_Y".format(bpm_name)] = bpm_data["{}_fit_Y".format(bpm_name)].apply(lambda x: np.abs(x[0][2]))

    return bpm_data


def bpm_plot(r, bpm_name, **kwargs):
    """Process a row of a BPM dataframe and returns a plot"""
    if kwargs.get("ax"):
        ax = kwargs.get("ax")
    else:
        ax = plt.figure().add_subplot(111)
    if kwargs.get("axis") == 'X' or kwargs.get("axis") == 'Y':
        axis = kwargs.get("axis")
    else:
        raise Exception("'axis' keyword argument must be 'X' or 'Y'.")

    x = r[f"{bpm_name}_data"]['channel']['data'][axis][:, 1]
    y = r[f"{bpm_name}_data"]['channel']['data'][axis][:, 2]

    ax.plot(x, y, 'b+:', label='data')
    ax.plot(x, gaussian(x, *r["{}_fit_{}".format(bpm_name, axis)][0]), 'ro:', label='fit')
