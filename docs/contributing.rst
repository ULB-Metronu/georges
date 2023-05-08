************
Contributing
************

Georges is an open-source Python library, and we encourage users to contribute to this repository.
The best way to contribute to the code is to make a fork of the repository and open a pull-request
to integrate the developments inside the master branch. To install all the packages, simply run ::

    poetry install --with dev,docs

Git Convention
##############

Principles
----------
* No one should commit to the master branch except the release manager.
* The master branch is the default people will get when cloning the repository and should only be a working tagged version. Therefore, the only commits in master are tagged versions.
* develop is our main development branch and features are made off of this branch.
* A new pull request must pass all tests before be merge inside the develop branch.
* Documentation should be added to the manual before merging the branch.

Naming
------
Feature branches should be descriptive and preferably named in the format “some-new-feature”.

Coding Style
############
We use `black <https://black.readthedocs.io/en/stable/>`_,
`isort <https://pycqa.github.io/isort/>`_ and `flake8 <https://flake8.pycqa.org/en/latest/>`_
to ensure consistent Python code. All the parameters are described in the `pyproject.toml` file.

Testing and metrics
###################
All the tests are available in the georges's test repository. We use `pytest <https://docs.pytest.org/en/7.2.x/>`_ as a tool for testing
and the results are then passed to `coverage` and uploaded to `sonar <https://www.sonarsource.com/products/sonarcloud/>`_
to compute metrics such as:

* Bugs
* Code smells
* Coverage

Documentation
#############

The documentation files are located in the docs repository and can be generated using sphinx. To create the HTML (or PDF) files,
go to the docs repository and use the `make` command::

    cd docs
    make html

The results are in the folder `build/html`.

Pre-commit file
###############

We provide a `.precommit.yml` file. This script simple can be executed before committing a file to automatically point out issues in
code such as missing semicolons, trailing whitespace, and debug statements.
`pre-commit <https://pre-commit.com>`_ can be installed using pip or brew and must be installed in the repository with::

    pre-commit install

And the script can be run using::

    pre-commit run --all-files -v

