# Installation

The installation process to get Georges up and running is relatively simple: the whole library is ready to be installed with `pip` and with `conda`.

However, as the :mod:`georges` submodules (in particular the :mod:`georges.manzoni` tracking code) are optimized to work with the [Intel distribution for Python](https://software.intel.com/en-us/articles/using-intel-distribution-for-python-with-anaconda), detailed installation instructions are provided to install :mod: georges with Conda. This is the recommended way for the installation.

The Georges' git repository is hosted on [Github](https://github.com/chernals/georges) and on [Gitlab](http://gitlab.sw.goiba.net/chernal/georges) (IBA internal).

## Obtaining the Georges source code with Git

Simply clone the repository:

    git clone https://github.com/chernals/georges

Or internally from IBA (replace _username_ with your SSO login):

    git clone git@gitlab.sw.goiba.net:username/georges.git

You can either stay on the bleeding-edge `master` branch or you can checkout a release tag:

    git checkout tags/2018.2

## Dependencies

A coherent set of dependencies is listed in the Conda `environment.txt` file as well as in the `setup.py` file (for `setuptools` and `pip`). In this way, installation using either `conda` or `pip` is possible. Note that the `pip` requirement file `requirements.txt` contains a single dot `.`, which refers to the dependency list provided in the `setup.py` file.

A typical user should not worry about those dependencies: they are properly managed either with `conda` (see next section) or with `pip` (see below).

## Installation with Anaconda and the Intel Python Distribution libraries

The installation procedure which follows creates a `conda` environment based on the Intel distribution for Python with all the necessary dependencies included and managed via `conda` itself (this ensures that all the dependencies are, if possible, based on the Intel channel and not managed with `pip`). Georges is then installed using `pip` from that `conda` environment. The dependencies are coherent and the `pip` installation of Georges will find all the dependencies listed in `setup.py` to be already installed in the `conda` environment.

1. Install [Conda](https://conda.io/docs/) for your operating system, follow the instructions for
  * [Linux](https://conda.io/docs/user-guide/install/linux.html)
  * [macOS](https://conda.io/docs/user-guide/install/macos.html)
  * [Windows](https://conda.io/docs/user-guide/install/windows.html)

2. Obtain a copy of the git repository (see previous section)

3. Create a dedicated `conda` environment (default name is `ipy3` for Intel Python 3)

        cd path/to/georges
        conda env create -f environment.yml

4. Activate the environment

        conda activate ipy3

5. Install Georges using `pip` from the `conda` environment

        # Typical installation
        pip install . 

Georges can be subsequently updated by running

    cd path/to/georges
    git pull origin master
    pip install --upgrade georges


## Installation with pip

In case Georges needs to be installed with the system Python or a Python installation that is not managed with `conda`, the following simple steps can be followed to use `pip` directly. All the dependencies are then managed with `pip`. Please note that this will typically not enable the Intel optimization, resulting in slower code execution for the :mod: manzoni module.

1. Obtain a copy of the git repository (see previous section)

2. Install Georges with `pip`:

        # Typical installation
        pip install . 

        # Install with pip in editable mode to get access to the modifications on the git repository
        pip install -e .

Georges can be subsequently updated by running

    cd path/to/georges
    git pull origin master
    pip install --upgrade georges


## Using Georges with Jupyter Notebook

Georges can be used with Jupyter notebooks. No special care is needed.

If you installed Georges within the `conda` environment, simply run (note that it is not advised to put all your notebook within the `git` structure):

    cd somewhere/good/for/notebooks
    jupyter notebook


## Georges distribution with Docker
__TODO__

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


