import numpy as np
from ..model import ManzoniModel
from .. import manzoni


class Cost:
    def __init__(self, model):
        if not isinstance(model, ManzoniModel):
            raise Exception("'model' must be an instance of ManzoniModel.")
        self._model = model


class CostTransmission(Cost):
    def __init__(self, model, **kwargs):
        super().__init__(model)

    def __call__(self, x):
        obs = manzoni.ElementByElementObserver(element_by_element=lambda b: [
            b.shape[0],
        ], elements=self._model.elements, track_start=None)
        manzoni.adjust_line(self._model.beamline, self._model.variables, x)
        o = manzoni.manzoni.track(self._model.beamline, self._model.beam, observer=obs)
        return (self._model.beam.shape[0] - o.data[0, 0, 0]) / self._model.beam.shape[0]


class CostOptics(Cost):
    def __init__(self, model, sigma_x={}, sigma_y={}, symmetry={}, alpha_x={}, alpha_y={}, transmission={}):
        super().__init__(model)
        self._parameters = dict()
        self._parameters['sigma_x'] = sigma_x
        self._parameters['sigma_y'] = sigma_y
        self._parameters['symmetry'] = symmetry
        self._parameters['alpha_x'] = alpha_x
        self._parameters['alpha_y'] = alpha_y
        self._parameters['transmission'] = transmission

    def __call__(self, x):
        obs = manzoni.ElementByElementObserver(element_by_element=lambda b: [
            np.std(b[:, 0]),
            np.std(b[:, 2]),
            np.cov(b[:, 0:2].T)[0, 1],
            np.cov(b[:, 2:4].T)[0, 1],
            b.shape[0],
        ], elements=self._model.elements)
        manzoni.adjust_line(self._model.beamline, self._model.variables, x)
        o = manzoni.manzoni.track(self._model.beamline, self._model.beam, observer=obs)
        p = self._parameters
        return (
            p['sigma_x'].get('weight', 0.0) * (o.data[0, 0, 0] - p['sigma_x'].get('value', 0.0)) ** 2
            + p['sigma_y'].get('weight', 0.0) * (o.data[0, 0, 1] - p['sigma_y'].get('value', 0.0)) ** 2
            + p['symmetry'].get('weight', 0.0) * (
                np.abs(o.data[0, 0, 0] - o.data[0, 0, 1]) / (o.data[0, 0, 0] + o.data[0, 0, 1]))
            + p['alpha_x'].get('weight', 0.0) * (o.data[0, 0, 2] - p['alpha_x'].get('value', 0.0)) ** 2
            + p['alpha_y'].get('weight', 0.0) * (o.data[0, 0, 3] - p['alpha_y'].get('value', 0.0)) ** 2
            + p['transmission'].get('weight', 0.0) * (self._model.beam.shape[0] - o.data[0, 0, 4]) /
            self._model.beam.shape[0]
        )
