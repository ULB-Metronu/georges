"""
TODO
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Optional, List
from georges_core.sequences import Sequence as _Sequence
from . import elements
from .core import track
from .elements.scatterers import MaterialElement
if TYPE_CHECKING:
    from georges_core import ureg as _ureg
    from .beam import Beam as _Beam
    from .observers import Observer as _Observer


class Input:
    def __init__(self, sequence: Optional[List[elements.ManzoniElement]] = None, beam: Optional[_Beam] = None):
        self._sequence = sequence
        self._beam = beam

    @property
    def sequence(self):
        return self._sequence

    @property
    def beam(self):
        return self._beam

    def freeze(self):
        """
        Freezes all elements in the input sequence.

        Returns:
            `self` to allow method chaining
        """
        for e in self._sequence:
            e.freeze()
        return self

    def unfreeze(self):
        """
        Unfreezes all elements in the input sequence.

        Returns:
            `self` to allow method chaining
        """
        for e in self._sequence:
            e.unfreeze()
        return self

    def track(self, beam: _Beam, observer: Optional[_Observer] = None, check_apertures: bool = True) -> _Observer:
        """

        Args:
            beam:
            observer:
            check_apertures:

        Returns:
            the `Observer` object containing the tracking results.
        """
        track(self, beam, observer, check_apertures)
        return observer

    def adjust_energy(self, input_energy: _ureg.Quantity):
        current_energy = input_energy
        for e in self.sequence:
            if isinstance(e, MaterialElement):
                e.KINETIC_ENERGY = current_energy
                current_energy = e.degraded_energy

    @classmethod
    def from_sequence(cls,
                      sequence: _Sequence,
                      ):
        """
        Creates a new `Input` from a generic sequence from `georges_core`.

        Args:
            sequence:

        Returns:

        """
        input_sequence = list()
        for name, element in sequence.df.iterrows():
            element_class = getattr(elements, element['CLASS'])
            parameters = set(list(element.index.values)).intersection(element_class.PARAMETERS.keys())
            input_sequence.append(
                element_class(name, **element[parameters])
            )
        return cls(sequence=input_sequence)
