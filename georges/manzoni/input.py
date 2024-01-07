"""
The file `input.py` contains the main class that needs instantiated to build a Manzoni input.
It requires at least a beamline Sequence and a Beam instance to allow the tracking. Once instantiated,
the Input class has a track function, which is just a wrapper on the track function of the file `core.py`,
and can be directly called to perform the tracking through the beamline. There are also a couple of
useful functions to freeze, unfreeze, set or get the parameters of a given element of the sequence.
Finally, the “adjust_energy” function is important to calculate all the initial kinetic energies of
all the scatterers and degraders of the beamline, starting from the initial kinetic energy of the beam.
It is necessary to correct the particles' propagation through these elements during the tracking.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Union

import numpy as _np
import pandas as _pd
from georges_core import Kinematics as _Kinematics
from georges_core import ureg as _ureg
from georges_core.sequences import BetaBlock as _BetaBlock
from georges_core.sequences import Sequence as _Sequence

from ..fermi import materials
from . import elements
from .beam import Beam as _Beam
from .core import track, twiss
from .elements import ManzoniElement
from .elements.scatterers import MaterialElement
from .integrators import Integrator, MadXIntegrator
from .observers import Observer as _Observer

MANZONI_FLAVOR = {"Sbend": "SBend", "Rbend": "RBend"}


class Input:
    def __init__(
        self,
        sequence: Optional[List[elements.ManzoniElement]] = None,
        beam: Optional[_Beam] = None,
        mapper: Dict[str, int] = None,
    ):
        self._sequence = sequence
        self._beam = beam
        self.mapper = mapper

    @property
    def sequence(self):  # pragma: no cover
        return self._sequence

    @property
    def beam(self):  # pragma: no cover
        return self._beam

    def to_df(self):
        """

        Returns: A pandas.DataFrame of the sequence

        """

        _ = list(map(lambda e: _pd.Series(e.attributes), self.sequence))
        df = _pd.concat(_, axis=1).T
        df["CLASS"] = list(map(lambda e: e.__class__.__name__, self.sequence))
        return df.set_index("NAME")

    @property
    def df(self):  # pragma: no cover
        return self.to_df()

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

    def track(
        self,
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
        track(self, beam, observers, check_apertures_exit=check_apertures)
        if observers is not None:
            if len(observers) == 1:
                return observers[0]
            else:
                return observers

    def twiss(
        self,
        kinematics: _Kinematics,
        reference_particle: _np.ndarray = None,
        offsets=None,
        twiss_parametrization: bool = True,
        twiss_init: _BetaBlock = None,
    ) -> _pd.DataFrame:
        """

        Args:
            kinematics:
            reference_particle:
            offsets:
            twiss_parametrization:
            twiss_init:

        Returns:
            The dataframe with the Twiss functions at each element.

        """
        return twiss(self, kinematics, reference_particle, offsets, twiss_parametrization, twiss_init)

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

    def set_integrator(self, integrator: Integrator = MadXIntegrator):
        for i in range(len(self.mapper)):
            self.sequence[i].integrator = integrator

    # TODO: use method __setitem__ instead ?
    def set_parameters(self, element: str, parameters: Dict):
        # unfreeze the element
        self.sequence[self.mapper[element]].unfreeze()
        for param in parameters.keys():
            self.sequence[self.mapper[element]].__setattr__(param, parameters[param])
        self.sequence[self.mapper[element]].freeze()

    def get_parameters(self, element: str, parameters: Optional[Union[List, str]] = None):
        if parameters is None:
            parameters = self.sequence[self.mapper[element]].attributes
            return dict(zip(parameters, list(map(self.sequence[self.mapper[element]].__getattr__, parameters))))
        if isinstance(parameters, List):
            return dict(zip(parameters, list(map(self.sequence[self.mapper[element]].__getattr__, parameters))))
        if parameters is not None and isinstance(parameters, str):
            return self.sequence[self.mapper[element]].__getattr__(parameters)

    @classmethod
    def insert_thin_element(
        cls,
        sequence: Input = None,
        position: int = 0,
        thin_element: ManzoniElement = None,
    ) -> Input:
        """
        Insert a thin element (e.g fringes) in the sequence.

        Args:
            sequence: the manzoni sequence
            position: index of the position to insert in the sequence
            thin_element: element to insert

        Returns:
            A new instance of manzoni.Input
        """

        new_sequence = sequence.sequence.copy()
        new_sequence.insert(position, thin_element)
        element_mapper = {k.NAME: v for v, k in enumerate(new_sequence)}
        return cls(sequence=new_sequence, mapper=element_mapper)

    @classmethod
    def from_sequence(
        cls,
        sequence: _Sequence,
        from_element: str = None,
        to_element: str = None,
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
        if "MATERIAL" in df_sequence.columns:
            idx = df_sequence[sequence.df["MATERIAL"].notnull()].index
            for ele in idx:
                if not isinstance(df_sequence.loc[ele, "MATERIAL"], materials.CompoundType):
                    df_sequence.loc[ele, "MATERIAL"] = getattr(materials, df_sequence.loc[ele, "MATERIAL"])

        for name, element in df_sequence.iterrows():
            element = element.replace({_np.nan: None})
            element_class = getattr(elements, MANZONI_FLAVOR.get(element["CLASS"], element["CLASS"]))
            parameters = list(set(list(element.index.values)).intersection(element_class.PARAMETERS.keys()))
            params = element[parameters].dropna()  # Remove the None from the parameters.
            input_sequence.append(
                element_class(name, **params),
            )
        element_mapper = {k: v for v, k in enumerate(list(df_sequence.index.values))}
        return cls(sequence=input_sequence, mapper=element_mapper)
