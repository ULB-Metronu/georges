from . import manzoni


class Model:
    def __init__(self, model=None, beam=None, beamline=None, context={}, variables=[], elements=[]):
        if model is not None and isinstance(model, Model):
            self._beam = model.beam
            self._beamline = model.beamline
            self._context = model.context
            self._variables = model.variables
            self._elements = model.elements
        else:
            self._beam = beam
            self._beamline = beamline
            self._context = context
            self._variables = variables
            self._elements = elements

    @property
    def beam(self):
        return self._beam

    @property
    def beamline(self):
        return self._beamline

    @property
    def context(self):
        return self._context

    @context.setter
    def context(self, c):
        self._context = c

    @property
    def variables(self):
        return self._variables

    @property
    def elements(self):
        return self._elements

    @property
    def gbeam(self):
        return self._beam

    @property
    def gbeamline(self):
        return self._beamline

    @property
    def gvariables(self):
        return self._variables

    @property
    def gelements(self):
        return self._elements


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
