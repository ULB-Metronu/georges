================
Model Comparison
================

Once a BDSIM model has been prepared from another model, it is of interest
to validate it to ensure the model has been prepared correctly.



Preparing Optics with BDSIM
---------------------------

The BDSIM model should be run with a 'core' beam distribution - ie typically
a Gaussian or Twiss Gaussian that will match the optics of the lattice. For
a physics study one might use a halo, but this is unsuitable for optics validation.

To compare, a BDSIM model is run with samplers attached to each element. This
records all of the particle coordinates at the end of each element. Once finished
a separate program ('rebdsim') is used to calculate moments and optical functions
from the distribution at each plane. This information can then be compared to
an anlytical description of the lattice such as that from MADX.

.. note:: It is important to open any apertures that are by design close to the beam
	  such as collimators. A non-Gaussian distribution will affect the calculation of
	  the optical parameters from the particle distribution.

Running BDSIM
*************

We recommend the following settings:

* Collimators are opened to at least 6 sigma of the beam distribution at their location.
* The `stopSecondaries` and `stopTracks` options are turned on to prevents secondaries being
  simulated and recorded.
* The physics list is set to "" - an empty string. This leaves only magnetic field tracking so
  that if a particle does hit the accelerator it will pass through without scattering.
* Simulate between 1000 and 50000 particles (events).

.. note:: This procedure is only suited to comparing linear optical functions. If sextupoles
	  or higher order magnets are present, these should be set to zero strength but must
	  remain in the lattice. The pybdsim.Convert.MadxTfs2Gmad converter for example provides
	  a boolean flag to convert the lattice with only linear optical components. The user
	  may of course proceed with non-linear magnetic fields included but it is only useful
	  to compare the sigma in each dimension to a similarly similuated distribution and not
	  the Twiss parameters.


Analysing Optical Data
**********************

The `rebdsim` tool can be used with an input `analysisConfig.txt` that specifies
`CalculateOpticalFunctions` to 1 or true in the header (see BDSIM manual). Or
the specially prepared optics tool `rebdsimOptics` can be used to achieve the
same outcome - we recommend this. In the terminal::

  $> rebdsimOptics myOutputFile.root optics.root

This may take a few minutes to process. This analyses the file from the BDSIM run
called 'myOutputFile.root' and produces another ROOT file called `optics.root` with
a different structure. This output file contains only optical data.

Comparing to MADX
-----------------

After preparing the optics from BDSIM, they may be compared to a MADX Tfs instance
with the following command in Python (for example)::

  >>> pybdsim.Compare.MadxVSBDSIM('twiss_v5.2fs', 'optics.root')

This will produce a series of plots comparing the orbit, beam size, and linear
optical functions.

The MADX twiss file (in tfs format) should contain all the possible columns in
the Twiss Module table. This can be prepared in a similar way as we would do
for converting to BDSIM GMAD syntax::

  select,flag=twiss, clear;
  twiss,sequence=SEQUENCENAME, file=twiss.tfs;


Comparing to MAD8
-----------------

Comparing to Transport
----------------------
