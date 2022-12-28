************
Installation
************

Georges is hosted on Github and can be downloaded using the following::

    git clone https://github.com/ULB-Metronu/georges.git

You can  stay on the bleeding-edge master branch or you can checkout
a release tag::

    git checkout tags/2023.1

The installation process to get George’s running is relatively simple: the whole library is ready to
be installed with `Poetry <https://python-poetry.org/>`_.

Dependencies
############

A coherent set of dependencies is listed in the `pyproject.toml` file (section tool.poetry.dependencies)
A typical user should not worry about those dependencies: they are correctly managed either with poetry

Installation using Poetry
#########################

Assuming you have Poetry and Python installed on your system, go to the location of the library and simply use
these commands::

    cd path/to/georges
    poetry install --without dev,docs

.. note::

    George’s uses python version >=3.8.1 and < 3.11

Georges can be subsequently updated by running the following::

    git pull origin master
    poetry update

Using Georges with Jupyter Lab
##############################

Georges can be used with Jupyter lab. No special care is needed,
and you can simply run (note that it is not advised to put all your
notebook within the git structure)::

    cd somewhere/good/for/notebooks
    jupyter-lab


Georges distribution with Docker
################################

A Docker image is made available to provide an easy access to a
complete Jupyter Lab + georges environment.

Use the *Dockerfile* to build the image::

    docker build

or, to register the image as well::

    docker build -t georges- -fDockerfile .

You can run a container with::

    docker run -it --rm --name georges -p 8899:8899 georges

then connect to http://127.0.0.1:8899 to access the Jupyter Lab interface
and type::

    import georges
