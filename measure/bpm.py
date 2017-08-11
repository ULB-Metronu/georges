import re
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


DATA_BPM = {
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