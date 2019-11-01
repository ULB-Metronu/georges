"""
TODO
"""
from typing import Optional, List
from georges_core.sequences import Sequence as _Sequence
from . import elements


class Input:
    def __init__(self, sequence: Optional[List[elements.ManzoniElement]] = None):
        self._sequence = sequence

    @property
    def sequence(self):
        return self._sequence

    def freeze(self):
        for e in self._sequence:
            e.freeze()
        return self

    def unfreeze(self):
        for e in self._sequence:
            e.unfreeze()
        return self

    @classmethod
    def from_sequence(cls,
                      sequence: _Sequence,
                      ):
        input_sequence = list()
        for name, element in sequence.df.iterrows():
            element_class = getattr(elements, element['CLASS'])
            parameters = set(list(element.index.values)).intersection(element_class.PARAMETERS.keys())
            input_sequence.append(
                element_class(name, **element[parameters])
            )
        return cls(sequence=input_sequence)
