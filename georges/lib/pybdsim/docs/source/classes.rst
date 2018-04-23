===============
Utility Classes
===============

Various classes are provided for the construction of BDSIM input blocks. Each class
can be instantiated and then used to prepare the gmad syntax using the Python
`str` or `repr` functions. These are used by the builder classes as well as the
converter functions.

Beam.Beam
---------

This beam class represents a beam definition in gmad syntax. The class has 'setter'
functions that are added dynamically based on the distribution type selected.::

  >>> b = pybdsim.Beam.Beam()
  >>> b.SetParicleType("proton")
  >>> b.SetDistributionType("reference")
  

Field
-----

This module allows BDSIM format field maps to be written and loaded. There are also
some plotting functions.  Please see :ref:`pybdsim-field-module` for more details.

Options.Options
---------------

This class provides the set of options for BDSIM. Please see
:ref:`pybdsim-options-module` for more details.

XSecBias.XSecBias
-----------------

This class provides the definition process biasing in BDSIM. Please see
:ref:`pybdsim-xsecbias-module` for more details.
