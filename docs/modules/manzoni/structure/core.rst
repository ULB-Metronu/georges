.. core:

Core
----

The beam propagation from one element to another in the beamline is implemented in the file `core.py`.
The `track` function contains the loop over the different elements, with optionally the check of
the apertures to select the particles that survive the tracking at the end of each element,
and the use of “observers” if defined by the user, to save the data during the tracking.
The file `core.py` also contains a Twiss function, which allows the calculation of the matrix elements
of all the elements along a given beamline for the Twiss functions calculation based on the 11 particles
method. The user must be aware that this function also needs the georges_core module to work properly,
as the Twiss computation is done in the end in the georges_core library, using the matrix elements
calculated in georges.