"""
TODO
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Optional, List, Union
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

    def track(self,
              beam: _Beam,
              observers: Union[List[_Observer], _Observer] = None,
              check_apertures: bool = True,

              ) -> Union[List[_Observer], _Observer]:
        """

        Args:
            beam:
            observers:
            check_apertures:

        Returns:
            the `Observer` object containing the tracking results.
        """
        if not isinstance(observers, list):
            observers = [observers]
        track(self, beam, observers, check_apertures)
        if observers is not None:
            if len(observers) == 1:
                return observers[0]
            else:
                return observers

    def adjust_energy(self, input_energy: _ureg.Quantity):
        current_energy = input_energy
        for e in self.sequence:
            if not isinstance(e, MaterialElement):
                break
            e.KINETIC_ENERGY = current_energy
            current_energy = e.degraded_energy

    def compute_efficiency(self, input_energy: _ureg.Quantity) -> float:
        self.adjust_energy(input_energy)
        efficiency = 1.0
        for e in self._sequence:
            if isinstance(e, MaterialElement):
                efficiency *= e.cache[5]
        return efficiency

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
