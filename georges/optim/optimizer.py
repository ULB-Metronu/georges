import scipy.optimize
from .. import manzoni
from ..model import Model, ManzoniModel


class Optimizer:
    LOCAL_METHODS = ['nelder-mead', 'powell']
    GLOBAL_METHODS = ['bassinhopping', 'diffevolution', 'ga']

    def __init__(self, model=None, cost=None, disp=True, method=None, **kwargs):
        if not isinstance(model, Model):
            raise Exception("'model' must be an instance of Model.")
        if isinstance(model, ManzoniModel):
            raise Exception("'model' cannot be an instance of ManzoniModel.")
        self._model = model
        self._manzoni_model = ManzoniModel(model)
        self._cost = cost(self._manzoni_model, **kwargs)
        self._disp = disp
        self._method = method
        self._kwargs = kwargs
        if self._method is None:
            raise Exception(
                f"Invalid method; should be one of {set([Optimizer.LOCAL_METHODS, Optimizer.GLOBAL_METHODS])}")

    @property
    def result(self):
        return self._result

    @property
    def tracking_result(self):
        return self._tracking_result

    @property
    def method(self):
        return self._method

    def _compute_results(self):
        self._tracking_result = manzoni.track(line=self._model.beamline,
                                              context=get_optimized_context(self.result.x),
                                              beam=self._model.beam)

    def _local_minimize(self, x0):
        self._result = scipy.optimize.minimize(self._cost,
                                               x0=x0,
                                               method=self.method,
                                               options={'disp': True},

                                               )
        if self._disp is not None:
            if self._disp:
                print(self._result)

    def _bassin_hopping_minimize(self, x0):
        self._result = scipy.optimize.basinhopping(self._cost,
                                                   x0=x0,
                                                   niter=200,
                                                   T=0.005,
                                                   stepsize=0.75,
                                                   disp=True,
                                                   minimizer_kwargs={
                                                       'method': 'SLSQP',
                                                       'bounds': [[-20, 0], [0, 20], [0, 20], [-20, 0], [0, 20]],
                                                   }
                                                   )
        if self._disp is not None:
            if self._disp:
                print(self._result)

    def run(self, x0):
        if self._method in Optimizer.LOCAL_METHODS:
            self._local_minimize(x0)
        elif self._method in Optimizer.GLOBAL_METHODS:
            if self._method is 'bassinhopping':
                self._bassin_hopping_minimize(x0)
        if self.result.success:
            self._compute_results()
        return self
