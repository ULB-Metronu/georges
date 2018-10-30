from typing import Callable
import numpy as np

_ = lambda x: x.copy()


class Observer:

    def __init__(self, turns: int=1, elements=[], func: Callable=_):
        self._data = np.empty(shape=(turns, max((len(elements), 1))), dtype=object)
        self._turns = turns
        self._elements = elements
        self._func = func

    @property
    def data(self):
        return self._data

    @property
    def turns(self) -> int:
        return self._turns

    def _observe(self, turn, element, beam):
        self._data[turn, element] = self._func(beam)
