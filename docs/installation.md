# Installation

The installation process to get Georges up and running is relatively simple: the whole library is ready to be installed with `pip`.

However, as the :mod:`georges` submodules (in particular the :mod:`georges.manzoni` tracking code) are optimized to work with the [Intel distribution for Python ](https://software.intel.com/en-us/articles/using-intel-distribution-for-python-with-anaconda), detailed installation instructions are provided to install :mod: georges with Conda. This is the recommended way for the installation.


## Installation with Anaconda and the Intel Python Distribution libraries
The installation procedure which follows creates a Conda environment with all the necessary dependencies included and managed via conda itself. Georges is then installed using `pip` from that `conda` environment. The dependencies are coherent and the `pip` installation of Georges will find all the dependencies listed in `setup.py` to be already installed in the `conda` environment.

1. Install [Conda](https://conda.io/docs/) for your operating system, follow the instructions for
  * [Linux](https://conda.io/docs/user-guide/install/linux.html)
  * [macOS](https://conda.io/docs/user-guide/install/macos.html)
  * [Windows](https://conda.io/docs/user-guide/install/windows.html)

## Dependencies

A coherent set of dependencies is listed in the Conda `environment.txt` file as well as in the `setup.py` file (for `setuptools` and `pip`). In this way, installation using either `conda` or `pip` is possible.

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
 
## Using Georges with Jupyter Notebook
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


download anaconda

conda update conda

