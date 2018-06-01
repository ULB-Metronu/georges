import numpy as np
import scipy.optimize
from .. import manzoni
from ..model import Model, ManzoniModel


class Optimizer:
    LOCAL_METHODS = ['nelder-mead', 'powell']
    GLOBAL_METHODS = ['bassinhopping', 'diffevolution', 'ga']

    def __init__(self, model=None, cost=None, disp=True, debug=False, **kwargs):
        if not isinstance(model, Model):
            raise Exception("'model' must be an instance of Model.")
        if isinstance(model, ManzoniModel):
            raise Exception("'model' cannot be an instance of ManzoniModel.")
        self._model = model
        self._manzoni_model = ManzoniModel(model)
        self._cost = cost(self._manzoni_model, **kwargs)
        self._disp = disp
        self._debug = debug
        self._method = None
        self._stored_x = []
        self._stored_f = []
        self._stored_accept = []

    @property
    def result(self):
        return self._result

    @property
    def tracking_result(self):
        return self._tracking_result

    @property
    def model(self):
        return self._model

    @property
    def manzoni_model(self):
        return self._manzoni_model

    @property
    def storage(self):
        return list(map(np.array, [self._stored_x, self._stored_f, self._stored_accept]))

    @property
    def method(self):
        return self._method

    def _get_optimized_context(self):
        print(self._model.context)
        for i, v in enumerate(self._model.variables):
            self._model.context[self._model.beamline.line.at[v[0], 'CIRCUIT']] = self._result.x[i]
        return self._model.context

    def _storage_callback(self, x, f, accept):
        self._stored_x.append(x)
        self._stored_f.append(f)
        self._stored_accept.append(accept)

    def _compute_results(self):
        self._tracking_result = manzoni.track(line=self._model.beamline,
                                              context=self._get_optimized_context(),
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

    def _bassin_hopping_minimize(self,
                                 x0=None,
                                 niter=10,
                                 T=1.0,
                                 stepsize=1.0,
                                 minimizer='SLSQP',
                                 callback=None,
                                 accept_test=None,
                                 ):



        # Default value for the initial guess
        if x0 is None:
            x0 = np.zeros(len(self._model.variables))

        # Call the minimization method from scipy
        self._result = scipy.optimize.basinhopping(
            self._cost,
            x0=x0,
            niter=niter,
            T=T,
            stepsize=stepsize,
            disp=self._debug,
            minimizer_kwargs={
                'method': minimizer,
                'bounds': [[-20, 0], [0, 20], [0, 20], [-20, 0], [0, 20]],
            },
            callback=callback,
            accept_test=accept_test,
        )

        # Display results
        if self._disp is not None:
            if self._disp:
                print(self._result)

    def run(self, x0, method=None, store=False, **kwargs):
        self._method = method
        if self._method not in set(Optimizer.LOCAL_METHODS + Optimizer.GLOBAL_METHODS):
            raise Exception(
                f"Invalid method; should be one of {set(Optimizer.LOCAL_METHODS + Optimizer.GLOBAL_METHODS)}")

        # Local optimization methods
        if self._method in Optimizer.LOCAL_METHODS:
            self._local_minimize(x0)
            if self._result['success']:
                self._compute_results()

        # Global optimization methods
        elif self._method in Optimizer.GLOBAL_METHODS:
            # Bassin hopping
            if self._method is 'bassinhopping':
                if kwargs.get('callback') is None and store:
                    kwargs['callback'] = self._storage_callback
                self._bassin_hopping_minimize(x0, **kwargs)
                if self._result.lowest_optimization_result.success:
                    self._compute_results()

            # Differential evolution
            if self._method is 'diffevolution':
                return

            # Genetic algorithm
            if self._method is 'ga':
                return
