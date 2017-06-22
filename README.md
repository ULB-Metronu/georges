# georges
Georges' the lemur opinionated particle accelerator modeling Python package.

<img src="https://github.com/chernals/georges/blob/master/georges.png" width="300" />

## Design
The aim of this library is to unify the description and computation of particle accelerator beamlines for different tools (MAD-X, PTC and GBeamline at this stage) in a unique Python library.

The library design strongly follows conventions and uses of the *Pandas* library: beamlines are naturally described as *Dataframes*. A functional approach is also one of the design goal of the library, something which fits well with *Pandas*.

Beamlines are loaded, converted (if necessary) and then transformed using functions split in packages (one package per beam physics code, *e.g.* MAD-X or G4Geamline). Those functional packages are supported by a series of *proxy classes* for each external computation code.
 
Support tools are also provided, notably a plotting library (entirely based on *Matplotlib*) which provides plotting capabilities for various optics computation (beam envelope, Twiss parameters, etc.).

## Usage




## Docker image
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