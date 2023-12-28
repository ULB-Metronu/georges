# Georges

[![ci](https://github.com/ULB-Metronu/georges/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/ULB-Metronu/georges/actions/workflows/master.yml)
[![documentation](https://github.com/ULB-Metronu/georges/actions/workflows/documentation.yml/badge.svg?branch=master)](https://github.com/ULB-Metronu/georges/actions/workflows/documentation.yml)
![Python](docs/_static/python_versions.svg)
![version](https://img.shields.io/badge/version-2023.2-blue)

[![Bugs](https://sonarcloud.io/api/project_badges/measure?project=ULB-Metronu_georges&metric=bugs)](https://sonarcloud.io/summary/overall?id=ULB-Metronu_georges)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=ULB-Metronu_georges&metric=coverage)](https://sonarcloud.io/summary/overall?id=ULB-Metronu_georges)
[![Reliability Rating](https://sonarcloud.io/api/project_badges/measure?project=ULB-Metronu_georges&metric=reliability_rating)](https://sonarcloud.io/summary/overall?id=ULB-Metronu_georges)

[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)
[![Gitter](https://badges.gitter.im/ULB-Metronu/georges.svg)](https://gitter.im/ULB-Metronu/georges?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)


Georges' the lemur opinionated particle accelerator modeling Python package.

<img src="docs/_static/georges.png" alt="drawing" width="300"/>

## Design
`Georges` provides a formalism for propagating many particles through magnetic elements while considering energy degradation by implementing the Fermi-Eyges technique. Beamlines are loaded and converted using [Georges-core](https://ulb-metronu.github.io/georges-core/index.html), where support tools are also provided, notably a plotting library (entirely based on Matplotlib and Plotly). It provides plotting capabilities for various optics computations (beam envelope, Twiss parameters, ...). 
Additionally, `Georges` includes a module for analyzing Bragg Peaks and estimating clinical properties such as R90 or lateral penumbra. It is also possible to determine the weight of each Bragg Peak in order to compute a Spread Out Bragg Peak (SOBP).

## Installation
`georges` is available from PyPI with pip:

    pip install georges

For a custom installation, please read the installation section in the [documentation](https://ulb-metronu.github.io/georges/installation.html).
