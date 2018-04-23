============
Installation
============


Requirements
------------

 * pytransport was developed for the Python 2.7 series.

pytransport depends on the following Python packages not included with Python:

 * matplotlib
 * numpy
 * scipy
 * pymadx
 * pybdsim

Installation
------------

A `setup.py` file required for a correct python installation is currently under development.

Currently, we recommend the user clones the source repository and exports the parent directory
to their PYTHONPATH environmental variable. This will allow Python to find pytransport.::

  pwd
  /Users/nevay/physics/reps
  git clone http://bitbucket.org/jairhul/pytransport
  ls
  > pytransport
  export PYTHONPATH=/Users/nevay/physics/reps

  python
  >>> import pytransport # no errors!
