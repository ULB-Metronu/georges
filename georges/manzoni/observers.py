"""
TODO
"""
from typing import Optional, List
import numpy as _np
import pandas as _pd
from numba import njit
from lmfit import Model, Parameters


class ObserverType(type):
    pass


class Observer(metaclass=ObserverType):
    def __init__(self, elements):
        self.elements = elements
        self.data = []
        self.headers = []

    def __call__(self, element: Optional[List[str]], b1, b2):
        if self.elements is None:
            return True
        if element.NAME not in self.elements:
            return False
        else:
            return True

    def to_df(self) -> _pd.DataFrame:
        return _pd.DataFrame(self.data, columns=self.headers).set_index('NAME')


class BeamObserver(Observer):
    def __init__(self, elements: Optional[List[str]] = None, with_input_beams: bool = False):
        super().__init__(elements)
        self._with_input_beams = with_input_beams
        self.headers = ('NAME', 'BEAM_IN', 'BEAM_OUT')

    def __call__(self, element, b1, b2):
        if super().__call__(element, b1, b2):
            self.data.append((element.NAME, _np.copy(b1) if self._with_input_beams else None, _np.copy(b2)))


class SuperObserver(Observer):
    def __init__(self, elements: Optional[List[str]] = None):
        super().__init__(elements)
        self.headers = ('NAME',
                        )

    def __call__(self, element, b1, b2):
        if super().__call__(element, b1, b2):
            self.data.append((element.NAME,))


class MeanObserver(Observer):
    def __init__(self, elements: Optional[List[str]] = None):
        super().__init__(elements)
        self.headers = ('NAME',
                        'BEAM_IN_X',
                        'BEAM_OUT_X',
                        'BEAM_IN_Y',
                        'BEAM_OUT_Y',
                        'BEAM_IN_XP',
                        'BEAM_OUT_XP',
                        'BEAM_IN_YP',
                        'BEAM_OUT_YP',
                        'BEAM_IN_DPP',
                        'BEAM_OUT_DPP',
                        )

    def __call__(self, element, b1, b2):
        if super().__call__(element, b1, b2):
            self.data.append((element.NAME,
                              b1[:, 0].mean(),
                              b2[:, 0].mean(),
                              b1[:, 2].mean(),
                              b2[:, 2].mean(),
                              b1[:, 1].mean(),
                              b2[:, 1].mean(),
                              b1[:, 3].mean(),
                              b2[:, 3].mean(),
                              b1[:, 4].mean(),
                              b2[:, 4].mean(),
                              ))


class SigmaObserver(Observer):
    def __init__(self, elements: Optional[List[str]] = None):
        super().__init__(elements)
        self.headers = ('NAME',
                        'BEAM_IN_X',
                        'BEAM_OUT_X',
                        'BEAM_IN_Y',
                        'BEAM_OUT_Y',
                        'BEAM_IN_XP',
                        'BEAM_OUT_XP',
                        'BEAM_IN_YP',
                        'BEAM_OUT_YP',
                        'BEAM_IN_DPP',
                        'BEAM_OUT_DPP',
                        )

    def __call__(self, element, b1, b2):
        if super().__call__(element, b1, b2):
            self.data.append((element.NAME,
                              b1[:, 0].std(),
                              b2[:, 0].std(),
                              b1[:, 2].std(),
                              b2[:, 2].std(),
                              b1[:, 1].std(),
                              b2[:, 1].std(),
                              b1[:, 3].std(),
                              b2[:, 3].std(),
                              b1[:, 4].std(),
                              b2[:, 4].std(),
                              ))


class LossesObserver(Observer):
    def __init__(self, elements: Optional[List[str]] = None):
        super().__init__(elements)
        self.headers = ('NAME',
                        'PARTICLES_IN',
                        'PARTICLES_OUT',
                        'TRANSMISSION',
                        'LOSSES',
                        )

    def __call__(self, element, b1, b2):
        if super().__call__(element, b1, b2):
            self.data.append((element.NAME,
                              b1.shape[0],
                              b2.shape[0],
                              100 * (b2.shape[0] / b1.shape[0]),
                              100 * (1 - b2.shape[0] / b1.shape[0]),
                              ))

    def compute_global_transmission(self, global_transmission: float = 1.0):
        for elem in self.data:
            global_transmission *= elem[3] / 100
        return global_transmission * 100

    def compute_global_losses(self):
        return 100 - self.compute_global_transmission()

    # def adjust_for_efficiency(self, efficiency: float = 1.0):
    # do something with self.data
    # self.data['LOSSES'] /= efficiency
    # self.data['TRANSMISSION'] *= effiency
    #   ...


class SymmetryObserver(Observer):
    def __init__(self):
        super().__init__()

        self.headers = ('NAME',
                        'SYM_IN',
                        'SYM_OUT',
                        )

    def __call__(self, element, b1, b2):
        self.data.append((element.NAME,
                          abs(b1[:, 0].std() - b1[:, 2].std()) / (b1[:, 0].std() + b1[:, 2].std()),
                          abs(b2[:, 0].std() - b2[:, 2].std()) / (b2[:, 0].std() + b2[:, 2].std()),
                          ))


class IbaBpmObserver(Observer):
    def __init__(self, elements: Optional[List[str]] = None):
        super().__init__(elements)
        self.headers = ('NAME',
                        'BEAM_OUT_X',
                        'BEAM_OUT_Y',
                        )

    @staticmethod
    def fit_bpm(distribution: _np.array):

        @njit
        def gaussian(x, a, mu, sigma):
            return a * _np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)) / (_np.sqrt(2 * _np.pi) * sigma)

        def fit_bpm(d, maxfev=1e7):
            bs = _np.array(
                [-31, -19.8, -15.8, -11.8, -7.8, -5.8, -3.8, -1.8, 0.0, 1.8,
                 3.8, 5.8, 7.8, 11.8, 15.8, 19.8, 31]) / 1000
            bsp = (bs[1:] + bs[:-1]) / 2
            w = 1.0 / (bs[1:] - bs[:-1])
            w[0] *= 0.7
            w[-1] *= 0.7
            hist = _np.histogram(d, bs)
            x = bsp
            y = w * hist[0]
            ar = _np.trapz(y / _np.sum(y) * len(y), x)
            mean = _np.mean(x * y / _np.sum(y) * len(y))
            rms = _np.std(x * y / _np.sum(y) * len(y))

            params = Parameters()  # Instanciate the Parameters class, then add variables as keywords
            params.add('a', value=ar, min=-1e6, max=1e6)
            params.add('mu', value=mean, min=mean-_np.abs(mean), max=mean+_np.abs(mean))
            params.add('sigma', value=rms, min=0.5*rms, max=2*rms)
            params.add('max_nfev', value=maxfev)

            gmodel = Model(gaussian)
            result = gmodel.fit(data=y,
                                x=x,
                                params=params)

            # Print the fit result
            print('bpm fit result:', result.best_values['a'], result.best_values['mu'], result.best_values['sigma'])
            return [result.best_values['sigma']]

        return fit_bpm(distribution)

    def __call__(self, element, b1, b2):
        if super().__call__(element, b1, b2):
            if element.CLASS == 'Marker':
                print(element.NAME)  # To identify the BPM whose data are being fitted
                self.data.append((element.NAME,
                                  self.fit_bpm(b2[:, 0])[0],
                                  self.fit_bpm(b2[:, 2])[0],
                                  ))
