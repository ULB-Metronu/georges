# Installation

:mod:`georges` submodules (in particular the :mod:`georges.manzoni` tracking code) are optimized to work with 

[Using Intel distribution for Python with Anaconda](https://software.intel.com/en-us/articles/using-intel-distribution-for-python-with-anaconda)


## Dependencies

Pandas >= 0.21 required.

## Installation with Anaconda


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

then connect to [http://localhost:8888](http://localhost:8888) to access the Jupyter Notebook interface.

The image includes a complete Anaconda Python3 environment with the most common packages. 
The latest *MAD-X* development release is available in */usr/local/bin/madx*.