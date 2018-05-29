import numpy as np


class GenericBounds:
    def __init__(self, xmax=None, xmin=None):
        self.xmax = None
        self.xmin = None
        if xmax is not None:
            self.xmax = np.array(xmax)
        if xmin is not None:
            self.xmin = np.array(xmin)

    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        if self.xmax is None and self.xmin is None:
            self.xmax = np.ones(x.shape[0])
            self.xmin = -np.ones(x.shape[0])
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin
