import numpy as np
from .transfer import transfer
from . import *


def convert_line(line):
    def class_conversion(e):
        if e['CLASS'] == 'DRIFT':
            e['CLASS_CODE'] = CLASS_CODE_DRIFT
        elif e['CLASS'] == 'QUADRUPOLE':
            e['CLASS_CODE'] = CLASS_CODE_QUADRUPOLE
        else:
            e['CLASS_CODE'] = 10
        return e

    def apertype_conversion(e):
        if e['APERTYPE'] == 'CIRCLE':
            e['APERTYPE_CODE'] = APERTYPE_CODE_CIRCLE
        elif e['APERTYPE'] == 'RECTANGLE':
            e['APERTYPE_CODE'] = APERTYPE_CODE_RECTANGLE
        return e
    line = line.apply(class_conversion, axis=1)
    line = line.apply(apertype_conversion, axis=1)
    return line


def aperture_check(b, e):
    # Circular aperture
    if e[3] == 1:
        s = (b[:, 0]**2 + b[:, 2]**2) < e[4]**2
    # Rectangular aperture
    elif e[3] == 2:
        s = (b[:, 0] < e[4]) & (b[:, 2] < e[4])
    # Unknown aperture type
    else:
        return b
    return np.compress(s, b, axis=0)


def track(l, b):
    line = l[['CLASS_CODE', 'L', 'K1', 'APERTYPE_CODE', 'APERTURE']].as_matrix()
    r = range(0, line.shape[0])
    beams = []
    for i in r:
        matrix = transfer[int(line[i, 0])](line[i])
        if matrix is not None:
            b = b.dot(matrix.T)
        b = aperture_check(b, line[i])
        beams.append(b)
    return beams
