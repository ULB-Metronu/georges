import numpy as np

_ = lambda x: x.copy()


class Observer:

    def __init__(self, turns=1, elements=[], func=_):
        self._start_data = None
        self._end_data = None
        self._data = np.empty(shape=(turns, max((len(elements), 1))), dtype=object)
        self._turns = turns
        self._elements = elements
        self._func = func

    @property
    def start_data(self):
        return self._start_data

    @property
    def end_data(self):
        return self._end_data

    @property
    def data(self):
        return self._data

    @property
    def func(self):
        return self._func

    @property
    def turns(self):
        return self._turns

    @property
    def elements(self):
        return self._elements

    def track_start(self, beam):
        self._start_data = self.func(beam)

    def observe(self, turn, element, beam):
        self._data[turn, element] = self.func(beam)
