============
Data Loading
============

Utilies to load BDSIM output data. This is intended for optical function plotting
and small scale data extraction - not general analysis of BDSIM output.


Loading ROOT Data
-----------------

The output optics in the ROOT file from `rebdsim` or `rebdsimOptics` may be loaded
with pybdsim providing the `root_numpy` package is available.::

  >>> d = pybdsim.Data.Load("optics.root")
