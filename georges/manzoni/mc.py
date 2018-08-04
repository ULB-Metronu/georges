import numpy as np
from .constants import *

try:
    import numpy.random_intel as nprandom
except ModuleNotFoundError:
    import numpy.random as nprandom


def scatterer(e, b, **kwargs):
    b[:,1] += nprandom.normal(0.0, e[INDEX_FE_A0], int(b.shape[0]))
    b[:,3] += nprandom.normal(0.0, e[INDEX_FE_A0], int(b.shape[0]))
    return b


def degrader(e, b, **kwargs):
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
    CLASS_CODES['DEGRADER']: degrader,
    CLASS_CODES['SCATTERER']: scatterer,
}
