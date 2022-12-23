.. maps:

Maps
----

The `maps` submodule contains the physics of the propagation of the particles through each of the
beamline elements described in the `elements` submodule. For each element, a first and second-order
type propagation is implemented, allowing the user to select the order of the tracking that is suitable
for his specific application. It is done via the selection of the integrator type when building the
beamline (see later). The `Transport-type`, `MadX-type`, and `Mad8-type` maps are implemented and
available. The user should be aware that the canonical variables of the particles are not the same
for the three different integrator types, so the definition of the beam must be done according to
the integrator to be consistent.