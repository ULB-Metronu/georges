import numpy as np

_ = lambda x: x.copy()


class Observer:

    def __init__(self, turns=1, elements=[], func=_):
        self._data = np.empty(shape=(turns, max((len(elements), 1))), dtype=object)
        self._turns = turns
        self._elements = elements
        self._func = func

    @property
    def data(self):
        return self._data

    def _observe(self, turn, element, beam):
        self._data[turn, element] = self._func(beam)