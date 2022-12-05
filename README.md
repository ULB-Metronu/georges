# Georges

[![develop](https://github.com/ULB-Metronu/georges/actions/workflows/develop.yml/badge.svg?branch=develop)](https://github.com/ULB-Metronu/georges/actions/workflows/develop.yml)
[![develop](https://github.com/rtesse/georges/actions/workflows/documentation.yml/badge.svg?branch=develop)](https://github.com/ULB-Metronu/georges/actions/workflows/documentation.yml)
![Python](docs/_static/python_versions.svg)
![version](https://img.shields.io/badge/version-2022.1-blue)

[![Bugs](https://sonarcloud.io/api/project_badges/measure?project=rtesse_georges&metric=bugs)](https://sonarcloud.io/summary/new_code?id=rtesse_georges)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=rtesse_georges&metric=coverage)](https://sonarcloud.io/summary/new_code?id=rtesse_georges)
[![Reliability Rating](https://sonarcloud.io/api/project_badges/measure?project=rtesse_georges&metric=reliability_rating)](https://sonarcloud.io/summary/new_code?id=rtesse_georges)

[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![codecov](https://codecov.io/gh/ULB-Metronu/georges/branch/develop/graph/badge.svg?token=IN1M54718K)](https://codecov.io/gh/ULB-Metronu/georges)
![version](https://img.shields.io/badge/version-2019.1-blue)
[![Python 3.8](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/release/python-380/)

Georges' the lemur opinionated particle accelerator modeling Python package. Also a thin wrapper over MAD-X/PTC, BDSim and G4Beamline.

<img src="https://raw.githubusercontent.com/ULB-Metronu/georges/legacy/docs/_static/georges.png" alt="drawing" width="300"/>

## Design
The aim of this library is to unify the description and computation of particle accelerator beamlines for different tools (MAD-X, PTC, BDSim and G4Beamline at this stage) in a unique Python library.

The library design strongly follows conventions and typical uses of the *Pandas* library: beamlines are naturally described as *Dataframes*. A functional approach is also one of the design goal of the library, something which fits well with *Pandas*. A clear separation between the beamline data, the beamline computations (*e.g.* Twiss) and the computation context is achieved. Beamline data are immutable, the context (*e.g.* external parameters such as the energy, momentum offset, initial beam, etc.) is immutable and always passed explicitely as a parameter and the computational facilities are implemented follwing a modular and chainable functional approach, unifying the above mentioned aspects.

Beamlines are loaded, converted (if necessary) and then transformed using functions split in packages (one package per beam physics code, *e.g.* MAD-X or G4Geamline). Those functional packages are supported by a series of *proxy classes* for each external computation code.
 
Support tools are also provided, notably a plotting library (entirely based on *Matplotlib*) which provides plotting capabilities for various optics computation (beam envelope, Twiss parameters, etc.).
