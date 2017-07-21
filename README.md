# georges
Georges' the lemur opinionated particle accelerator modeling Python package. Also a thin wrapper over MAD-X/PTC, BDSim and G4Beamline.

<img src="https://github.com/chernals/georges/blob/master/georges.png" width="300" />

## Design
The aim of this library is to unify the description and computation of particle accelerator beamlines for different tools (MAD-X, PTC, BDSim and G4Beamline at this stage) in a unique Python library.

The library design strongly follows conventions and typical uses of the *Pandas* library: beamlines are naturally described as *Dataframes*. A functional approach is also one of the design goal of the library, something which fits well with *Pandas*. A clear separation between the beamline data, the beamline computations (*e.g.* Twiss) and the computation context is achieved. Beamline data are immutable, the context (*e.g.* external parameters such as the energy, momentum offset, initial beam, etc.) is immutable and always passed explicitely as a parameter and the computational facilities are implemented follwing a modular and chainable functional approach, unifying the above mentioned aspects.

Beamlines are loaded, converted (if necessary) and then transformed using functions split in packages (one package per beam physics code, *e.g.* MAD-X or G4Geamline). Those functional packages are supported by a series of *proxy classes* for each external computation code.
 
Support tools are also provided, notably a plotting library (entirely based on *Matplotlib*) which provides plotting capabilities for various optics computation (beam envelope, Twiss parameters, etc.).

## Usage
No attempt is made to support python versions earlier than CPython 3.5. Jython and alternative implementation have not been tested.

### Physics module ###
`georges.physics` is a simple standalone helper module providing various relativistic conversion functions as well as range computations for protons in water.

It is automatically importe when the main `georges` module is imported.

Pull requests are encouraged, as the module growths it will be necessary to split it in autonomous pieces.

```
import georges
georges.physics.energy_to_beta(0.1)
georges.physics.range_to_energy(4.0)
```


## Docker image (experimental and unstable)
A Docker image is made available to provide an easy access to a complete Jupyter Notebook + madx + georges environment.
 
Use  the *Dockerfile* to build the image:
 
```
docker build .
```

or, to register the image as well:

```
docker build -t username/georges .
```

You can run a container with

```
docker run -it username/georges
```

then connect to [http://localhost:8888](http://localhost:8888 "Jupyter Notebook") to access the Jupyter Notebook interface.

The image includes a complete Anaconda Python3 environment with the most common packages. 
The latest *MAD-X* development release is available in */usr/local/bin/madx*.