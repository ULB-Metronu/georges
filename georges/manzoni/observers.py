"""
TODO
"""
import numpy as _np
import pandas as _pd


class ObserverType(type):
    pass


class Observer(metaclass=ObserverType):
    def __init__(self):
        self.data = []
        self.headers = []

    def __call__(self, element, b1, b2):
        pass

    def to_df(self) -> _pd.DataFrame:
        return _pd.DataFrame(self.data, columns=self.headers)


class BeamObserver(Observer):
    def __init__(self):
        super().__init__()
        self.headers = ('LABEL1', 'BEAM_IN', 'BEAM_OUT')

    def __call__(self, element, b1, b2):
        self.data.append((element.LABEL1, b1, b2))


class MeanObserver(Observer):
    pass


class SigmaObserver(Observer):
    def __init__(self):
        super().__init__()
        self.headers = ('LABEL1', 'BEAM_IN_X', 'BEAM_OUT')

    def __call__(self, element, b1, b2):
        self.data.append((element.LABEL1, b1[0].std(), b2[0].std()))


class LossesObserver(Observer):
    pass
