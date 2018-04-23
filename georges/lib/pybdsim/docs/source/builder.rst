===============
Building Models
===============


pybdsim provides a series of classes that allows a BDSIM model to be built
programmatically in Python and finally written out to BDSIM input syntax
('gmad').


Creating A Model
----------------

The :code:`Machine` class provides the functionality to create a BDSIM model.
This would instantiated and a sequence is defined by adding accelerator
elements in order to that instance by calling functions such as
:code:`AddDipole()`. Extra information can then be associated with that
Machine instance and finally, it can be written out to a series of gmad
files as input to BDSIM. For example::

  >>> a = pybdsim.Builder.Machine()
  >>> a.AddDrift()


The arguments can generally be found by using a question mark on a function.::

  >>> a.AddDrift?
  Signature: a.AddDrift(name='dr', length=0.1, **kwargs)
  Docstring: Add a drift to the beam line
  File:      ~/physics/reps/pybdsim/Builder.py
  Type:      instancemethod

Adding Options
--------------

No options are required to run the most basic BDSIM model. However, it is often
advantageous to specify at leat a few options such as the physics list and default
aperture. To add options programmatically, there is an options class. This is
instantiated and then 'setter' methods are used to set values of parameters.
This options instance can then be assocated with a machine instance. For example::

  >>> o = pybdsim.Options.Options()
  >>> o.SetPhysicsList('em hadronic decay muon hadronic_elastic')

  >>> a = pybdsim.Builder.Machine()
  >>> a.AddOptions(o)


The possible options can be seen by using tab complete in ipython::

  >>> a.Set<tab>

.. note:: Only the most common options are currently implement. Please
	  see :ref:`feature-request` to request others.


Adding a Beam
-------------

A beam definition that specifies at least the particle type and total energy
is required to run a BDSIM model. The machine class will provide a default
such that the model will run 'out of the box', but is of course of interest
to specify these options. To add a beam definition, there is a beam class.
This is instantiated and then 'setter' methods are used to set values of
parameters. this beam instance can then be associated with a machine instance.
For example::

  >>> b = pybdsim.Beam.Beam()
  >>> b.SetDistributionType('reference')
  >>> b.SetEnergy(25, 'GeV')
  >>> b.SetParticleType('proton')

  >>> a = pybdsim.Builder.Machine()
  >>> a.AddBeam(b)

.. note:: More setter functions will dynamically appear based on the distribution
	  type set.

Writing a Machine
-----------------

Once completed, a machine can be written out to gmad files to be used as input
for BDSIM. This is done as follows::

  >>> a = pybdsim.Builder.Machine()
  >>> a.Write('outputfilename')

Units
-----

The user may supply units as strings that will be written to the gmad syntax
as a Python tuple. For example::

  >>> a = pybdsim.Builder.Machine()
  >>> a.AddDrift('d1', (3.2, 'm'))

This will result in the following gmad syntax::

  >>> print a[0]
  d1: drift, l=3.2*m;

.. note:: There is no checking on the string supplied, so it is the users
	  responsibility to supply a valid unit string that BDSIM will accept.

  
kwargs - Flexibility
--------------------

'kwargs' are optional keyword arguments in Python. This allows the user to
supply arbitrary options to a function that can be instpected inside the
function as a dictionary. BDSIM gmad syntax to define an element generally
follows the pattern::

  name : type, parameter1=value, parameter2=value;

Many parameters can be added and this syntax is regularly extended. It would
therefore be impractical to have every function with all the possible arguments.
To solve this problem, the :code:`**kwargs` argument allows the user to
specify any option that will be passed along and written to file in the element
definition as 'key=value'. For example::

  >>> a = pybdsim.Builder.Machine()
  >>> a.AddDrift('drift321', 3.2, aper1=5, aper2=4.5, apertureType="rectangular")
  
This will result in the following gmad syntax being written::

  >>> print a[0]
  drift321: drift, apertureType="rectangular", aper2=4.5, aper1=5, l=3.2;

Anywhere you see a function with the last argument as :code:`**kwargs`, this
feature can be used.

The arguments included in the function signatures are the minimum arguments
required for functionality.
