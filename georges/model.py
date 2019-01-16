from typing import Optional, Dict, List
import numpy as _np
from . import manzoni
from .beam import Beam
from .beamline import Beamline


class Model:
    def __init__(self, model=None, beam=None, beamline=None, context=None, variables=None, elements=None):
        if model is not None and isinstance(model, Model):
            self._beam = model.beam
            self._beamline = model.beamline
            self._context = model.context
            self._variables = model.variables
            self._elements = model.elements
        else:
            self._beam = beam
            self._beamline = beamline
            self._context = context or {}
            self._variables = variables or []
            self._elements = elements or []

    @property
    def beam(self) -> Beam:
        return self._beam

    @property
    def beamline(self) -> Beamline:
        return self._beamline

    @property
    def context(self) -> Dict:
        return self._context

    @context.setter
    def context(self, c):
        self._context = c

    @property
    def variables(self) -> List:
        return self._variables

    @property
    def elements(self):
        return self._elements

    def track(self, **kwargs):
        return manzoni.track(self, **kwargs)


class ManzoniModel(Model):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._manzoni_beamline = None
        self._manzoni_beamline_numpy = None
        self._manzoni_beam = None
        self._manzoni_variables = None
        self._manzoni_elements = None

    @property
    def beam(self):
        if self._manzoni_beam is None:
            self._manzoni_beam = self._beam.distribution.values
        return self._manzoni_beam

    def get_beamline(self, to_numpy=True):
        if to_numpy:
            if self._manzoni_beamline_numpy is None:
                self._manzoni_beamline_numpy = manzoni.convert_line(
                    self._beamline.line,
                    context=self.context,
                    to_numpy=to_numpy
                )
            return self._manzoni_beamline_numpy
        else:
            if self._manzoni_beamline is None:
                self._manzoni_beamline = manzoni.convert_line(
                    self._beamline.line,
                    context=self.context,
                    to_numpy=to_numpy
                )
            return self._manzoni_beamline

    beamline = property(get_beamline)

    @property
    def variables(self):
        if self._manzoni_variables is None:
            self._manzoni_variables = manzoni.transform_variables(self.get_beamline(to_numpy=False), self._variables)
        return self._manzoni_variables

    @property
    def elements(self):
        if self._manzoni_elements is None:
            self._manzoni_elements = manzoni.transform_elements(self.get_beamline(to_numpy=False), self._elements)
        return self._manzoni_elements

    def adjust_beamline(self, values):
        manzoni.adjust_line(self.beamline, self.variables, _np.array(values))

    def adjust_beam(self, beam):
        try:
            self._manzoni_beam = beam.distribution.values
        except AttributeError:
            self._manzoni_beam = beam

    def track(self, observer: Optional[manzoni.Observer] = None, **kwargs):
        observer = observer or manzoni.Observer(elements=self.elements)
        return manzoni.manzoni.track(self.beamline, self.beam, observer, **kwargs)
