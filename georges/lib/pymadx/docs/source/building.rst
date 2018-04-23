=================
Model Preparation
=================

pymadx contains a series of classes that can be used to programmatically
construct a MADX model. The main class is `pymadx.Builder.Machine`.::

  a = pymadx.Builder.Machine()
  a.AddDrift('drift1',1.3)
  a.AddQuadrupole('qf1',0.2,1.3454)
  a.Write('lattice1')

The functions available are documented in :ref:`pymadx-builder`, but can
also easily be found with the built in documentation::

  a = pymadx.Builder.Machine()
  a <tab>

to see the list of available functions. Each has a short description
and signature that can be viewed with a question mark.::

  a = pymadx.Builder.Machine()
  a.AddQuadrupole?
  Signature: a.AddQuadrupole(name='qd', length=0.1, k1=0.0, **kwargs)
  Docstring: <no docstring>
  File:      ~/physics/reps/pymadx/pymadx/Builder.py
  Type:      instancemethod


Aside from the lattice elements available, a `pymadx.Beam.Beam` instance
can be associated with the machine.::

  b = pymadx.Beam.Beam()
  a.AddBeam(b)
