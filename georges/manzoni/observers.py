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
    def __init__(self, with_input_beams: bool = False):
        super().__init__()
        self._with_input_beams = with_input_beams
        self.headers = ('LABEL1', 'BEAM_IN', 'BEAM_OUT')

    def __call__(self, element, b1, b2):
        self.data.append((element.LABEL1, b1.copy() if self._with_input_beams else None, b2.copy()))


class MeanObserver(Observer):
    def __init__(self):
        super().__init__()
        self.headers = ('LABEL1',
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
        self.data.append((element.LABEL1,
                          b1[0].mean(),
                          b2[0].mean(),
                          b1[2].mean(),
                          b2[2].mean(),
                          b1[1].mean(),
                          b2[1].mean(),
                          b1[3].mean(),
                          b2[3].mean(),
                          b1[4].mean(),
                          b2[4].mean(),
                          ))


class SigmaObserver(Observer):
    def __init__(self):
        super().__init__()
        self.headers = ('LABEL1',
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
        self.data.append((element.LABEL1,
                          b1[0].std(),
                          b2[0].std(),
                          b1[2].std(),
                          b2[2].std(),
                          b1[1].std(),
                          b2[1].std(),
                          b1[3].std(),
                          b2[3].std(),
                          b1[4].std(),
                          b2[4].std(),
                          ))


class LossesObserver(Observer):
    pass
