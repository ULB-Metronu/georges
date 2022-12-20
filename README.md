# Georges

[![ci](https://github.com/rtesse/georges/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/rtesse/georges/actions/workflows/develop.yml)
[![documentation](https://github.com/rtesse/georges/actions/workflows/documentation.yml/badge.svg?branch=develop)](https://github.com/ULB-Metronu/georges/actions/workflows/documentation.yml)
![Python](docs/_static/python_versions.svg)
![version](https://img.shields.io/badge/version-2022.1-blue)

[![Bugs](https://sonarcloud.io/api/project_badges/measure?project=rtesse_georges&metric=bugs)](https://sonarcloud.io/summary/new_code?id=rtesse_georges)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=rtesse_georges&metric=coverage)](https://sonarcloud.io/summary/new_code?id=rtesse_georges)
[![Reliability Rating](https://sonarcloud.io/api/project_badges/measure?project=rtesse_georges&metric=reliability_rating)](https://sonarcloud.io/summary/new_code?id=rtesse_georges)

[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)
[![Gitter](https://badges.gitter.im/ULB-Metronu/georges.svg)](https://gitter.im/ULB-Metronu/georges?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)


Georges' the lemur opinionated particle accelerator modeling Python package.

<img src="https://raw.githubusercontent.com/ULB-Metronu/georges/legacy/docs/_static/georges.png" alt="drawing" width="300"/>

## Design
The aim of this library is to unify the description and computation of particle accelerator beamlines for different tools (MAD-X, PTC, BDSim and G4Beamline at this stage) in a unique Python library.

* Fermi 
* Manzoni

Beamlines are loaded using *georges-core*. 
Support tools are also provided, notably a plotting library (entirely based on *Matplotlib*) which provides plotting capabilities for various optics computation (beam envelope, Twiss parameters, etc.).
