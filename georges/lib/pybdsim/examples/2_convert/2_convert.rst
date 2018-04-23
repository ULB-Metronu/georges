Convert
=======

1_madxtfs2gmad.py
-----------------

A sample `.tfs` file is provided that represents the first ~300 m of the
LHC beam 1, 3.5 TeV lattice, containing around 160 mixed elements. This
example loads the compressed (`.tar.gz`) tfs file and converts it to
BDSIM gmad syntax.

How to run::

  ./1_madxtfs2gmad.py
  bdsim --file=lhcb1.gmad

The optics and lattice from madx:

.. figure:: twiss35tevb1_short.pdf
	    :width: 70%

The converted lattice as visualised in BDSIM.
  
.. figure:: 1_madxtfs2gmad.png
	    :width: 70%
  
