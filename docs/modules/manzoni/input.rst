Input
#####

The file `input.py` contains the main class that needs instantiated to build a Manzoni input.
It requires at least a beamline Sequence and a Beam instance to allow the tracking. Once instantiated,
the Input class has a track function, which is just a wrapper on the track function of the file `core.py`,
and can be directly called to perform the tracking through the beamline. There are also a couple of
useful functions to freeze, unfreeze, set or get the parameters of a given element of the sequence.
Finally, the “adjust_energy” function is important to calculate all the initial kinetic energies of
all the scatterers and degraders of the beamline, starting from the initial kinetic energy of the beam.
It is necessary to correct the particles' propagation through these elements during the tracking.