"""
TODO
"""
import numpy as _np
import pandas as _pd


class ObserverType(type):
    pass


class Observer(metaclass=ObserverType):
    def __init__(self):
        self.data = {}

    def __call__(self, element, b1, b2):
        pass

    def to_df(self) -> _pd.DataFrame:
        return _pd.DataFrame(self.data)


class BeamObserver(Observer):
    def __call__(self, element, b1, b2):
        self.data[element.LABEL1] = (b1, b2)


class SigmaObserver(Observer):
    def __call__(self, element, b1, b2):
        self.data[element.LABEL1] = _np.array()
