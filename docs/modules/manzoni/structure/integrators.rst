.. integrators:

Integrators
-----------

The file `integrators.py` contains the definition of all the integrator types,
with a propagate function that is a wrapper to the `propagate` function of each of the
different elements in the beamline. Each integrator has a dictionary with the different elements
it can be used with, so that if this integrator is selected, the propagate function of the element
will be called via the propagate function of the integrator to allow the tracking.