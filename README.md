# georges
Georges' the lemur opinionated particle accelerator modeling Python package. Also a thin wrapper over MAD-X/PTC, BDSim and G4Beamline.

<img src="https://github.com/chernals/georges/blob/master/docs/_static/georges.png" width="300" />

## Design
The aim of this library is to unify the description and computation of particle accelerator beamlines for different tools (MAD-X, PTC, BDSim and G4Beamline at this stage) in a unique Python library.

The library design strongly follows conventions and typical uses of the *Pandas* library: beamlines are naturally described as *Dataframes*. A functional approach is also one of the design goal of the library, something which fits well with *Pandas*. A clear separation between the beamline data, the beamline computations (*e.g.* Twiss) and the computation context is achieved. Beamline data are immutable, the context (*e.g.* external parameters such as the energy, momentum offset, initial beam, etc.) is immutable and always passed explicitely as a parameter and the computational facilities are implemented follwing a modular and chainable functional approach, unifying the above mentioned aspects.

Beamlines are loaded, converted (if necessary) and then transformed using functions split in packages (one package per beam physics code, *e.g.* MAD-X or G4Geamline). Those functional packages are supported by a series of *proxy classes* for each external computation code.
 
Support tools are also provided, notably a plotting library (entirely based on *Matplotlib*) which provides plotting capabilities for various optics computation (beam envelope, Twiss parameters, etc.).
