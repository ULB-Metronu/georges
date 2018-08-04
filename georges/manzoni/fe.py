import numpy as np
from .constants import *

try:
    import numpy.random_intel as nprandom
except ModuleNotFoundError:
    import numpy.random as nprandom


def propagation_scatterer(e, b, **kwargs):
    b[X, X] += e[INDEX_FE_A2] + b[X, PX] * e[INDEX_LENGTH] + b[1, 1] * e[INDEX_LENGTH]**2
    b[X, PX] += e[INDEX_FE_A1] + b[PX, PX] * e[INDEX_LENGTH]
    b[PX, X] = b[X, PX]
    b[PX, PX] += e[INDEX_FE_A0]
    
    b[Y, Y] += e[INDEX_FE_A2] + b[Y, PY] * e[INDEX_LENGTH] + b[1, 1] * e[INDEX_LENGTH]**2
    b[Y, PY] += e[INDEX_FE_A1] + b[PY, PY] * e[INDEX_LENGTH]
    b[PY, Y] = b[Y, PY]
    b[PY, PY] += e[INDEX_FE_A0]

    return b


def propagation_degrader(e, b, **kwargs):
    b[PX, PX] += e[INDEX_FE_A0]
    b[PY, PY] += e[INDEX_FE_A0]
    return b


def mc_scatterer(e, b, **kwargs):
    b[:,1] += nprandom.normal(0.0, e[INDEX_FE_A0], int(b.shape[0]))
    b[:,3] += nprandom.normal(0.0, e[INDEX_FE_A0], int(b.shape[0]))
    return b


def mc_degrader(e, b, **kwargs):
    # Remove particles
    # if e[INDEX_FE_LOSS] != 0:
    #     idx = np.random.randint(b.shape[0], size=int((e[INDEX_FE_LOSS]) * b.shape[0]))
    #     b = b[idx, :]
    b += nprandom.multivariate_normal(
        [0.0, 0.0, 0.0, 0.0, 0.0],
        np.array(
            [
                [e[INDEX_FE_A2], e[INDEX_FE_A1], 0, 0, 0],
                [e[INDEX_FE_A1], e[INDEX_FE_A0], 0, 0, 0],
                [0, 0, e[INDEX_FE_A2], e[INDEX_FE_A1], 0],
                [0, 0, e[INDEX_FE_A1], e[INDEX_FE_A0], 0],
                [0, 0, 0, 0, e[INDEX_FE_DPP]]
            ]
        ),
        int(b.shape[0]))
    return b


mc = {
    CLASS_CODES['DEGRADER']: mc_degrader,
    CLASS_CODES['SCATTERER']: mc_scatterer,
}

fe = {
    CLASS_CODES['DEGRADER']: propagation_degrader,
    CLASS_CODES['SCATTERER']: propagation_scatterer,
}
