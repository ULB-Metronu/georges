# Installation

:mod:`georges` submodules (in particular the :mod:`georges.manzoni` tracking code) are optimized to work with 

[Using Intel distribution for Python with Anaconda](https://software.intel.com/en-us/articles/using-intel-distribution-for-python-with-anaconda)


## Python installation with Anaconda and the Intel libraries

### Dependencies

Pandas >= 0.21 required.

## Installation
The easiest way to use `georges` is to clone it directly from Github:

    cd ~/reps
    git clone https://github.com/chernals/georges
    cd georges
    git status
    
assuming you have a `reps` directory in your home directory.

From there it is easy to pull the latest updates:

    cd ~/reps/georges
    git pull
   
To be able to import the library you will need to have it in your python path. One way is to add it to the environment variable `$PYTHONPATH`. For example:

    export PYTHONPATH=$PYTHONPATH:~/reps/
echo $PYTHONPATH
    
This will add all the libraries in your `~/reps` directory to `$PYTHONPATH` and will make them available for import.

## Using and importing Georges in Python
You can access the library by simply importing it:

    import georges
   
This will include only the core components of Georges. The different Georges' modules must be imported separately, depending on your needs:

    import georges.madx
    import georges.bdsim
    import georges.manzoni
    import georges.plotting
    
See the examples below for a typical use case.

    import georges
    from georges.plotting import *
    import georges.manzoni
 
## Using Geoges with Jupyter Notebook
xxxx  


## Using the complete Georges tools from Docker image
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