"""
TODO
"""
from __future__ import annotations
from typing import Optional, List, Union, Dict
from georges_core.sequences import Sequence as _Sequence
from . import elements
from .core import track
from .elements.scatterers import MaterialElement
from ..fermi import materials

from georges_core import ureg as _ureg
from .beam import Beam as _Beam
from .observers import Observer as _Observer

MANZONI_FLAVOR = {"Sbend": "SBend"}


class Input:
    def __init__(self, sequence: Optional[List[elements.ManzoniElement]] = None,
                 beam: Optional[_Beam] = None,
                 mapper: Dict[str, int] = None):
        self._sequence = sequence
        self._beam = beam
        self._mapper = mapper

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
            if isinstance(e, MaterialElement):
                e.KINETIC_ENERGY = current_energy
                current_energy = e.degraded_energy

    def compute_efficiency(self, input_energy: _ureg.Quantity) -> float:
        self.adjust_energy(input_energy)
        efficiency = 1.0
        for e in self._sequence:
            if isinstance(e, MaterialElement):
                efficiency *= e.cache[5]
        return efficiency

    # TODO: use method __setitem__ instead ?
    def set_parameters(self, element: str, parameters: Dict):
        # unfreeze the element
        self.sequence[self._mapper[element]].unfreeze()
        for param in parameters.keys():
            self.sequence[self._mapper[element]].__setattr__(param, parameters[param])
        self.sequence[self._mapper[element]].freeze()

    def get_parameters(self, element: str, parameters: Optional[List] = None):
        if parameters is None:
            parameters = self.sequence[self._mapper[element]].attributes
        return dict(zip(parameters, list(map(self.sequence[self._mapper[element]].__getattr__, parameters))))

    @classmethod
    def from_sequence(cls,
                      sequence: _Sequence,
                      from_element: str = None,
                      to_element: str = None
                      ):
        """
        Creates a new `Input` from a generic sequence from `georges_core`.

        Args:
            sequence:
            from_element:
            to_element:
        Returns:

        """
        input_sequence = list()
        df_sequence = sequence.df.loc[from_element:to_element]
        if 'MATERIAL' in df_sequence.columns:
            idx = df_sequence[sequence.df['MATERIAL'].notnull()].index
            for ele in idx:
                df_sequence.loc[ele, "MATERIAL"] = getattr(materials, df_sequence.loc[ele, "MATERIAL"])

        for name, element in df_sequence.iterrows():
            element_class = getattr(elements, element['CLASS'])
            parameters = list(set(list(element.index.values)).intersection(element_class.PARAMETERS.keys()))
            input_sequence.append(
                element_class(name, **element[parameters])
            )
        element_mapper = {k: v for v, k in enumerate(list(df_sequence.index.values))}
        return cls(sequence=input_sequence, mapper=element_mapper)
