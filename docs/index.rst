:tocdepth: 2

###################################
Welcome to Georges's documentation!
###################################

|Actions Status| |Documentation Status| |Python version| |version| |Bugs| |Coverage| |Reliability|
|License| |Code Style| |Gitter|


:code:`Georges` provides a formalism for propagating many particles through magnetic elements while considering energy degradation by implementing the Fermi-Eyges technique. Beamlines are loaded and converted using `Georges-core <https://ulb-metronu.github.io/georges-core/index.html>`_,  where support tools are also provided, notably a plotting library (entirely based on Matplotlib and Plotly). It provides plotting capabilities for various optics computations (beam envelope, Twiss parameters, ...). 
Additionally, :code:`Georges` includes a module for analyzing Bragg Peaks and estimating clinical properties such as R90 or lateral penumbra. It is also possible to determine the weight of each Bragg Peak in order to compute a Spread Out Bragg Peak (SOBP).

***********************
Georges's documentation
***********************

The documentation is part of the Georges-core repository itself and is made available *via* `Github Pages <https://pages.github.com>`_ . It is hosted at `ulb-metronu.github.io/georges/ <https://ulb-metronu.github.io/georges/>`_

We value your contributions and you can follow the instructions in :doc:`Contributing <contributing>`.
Finally, if you’re having problems, please do let us know at our :doc:`Support <support>` page.

..  toctree::
    :maxdepth: 2
    :titlesonly:
    :caption: Developers
    :glob:

    authors
    support
    contributing
    changelog

..  toctree::
    :maxdepth: 2
    :titlesonly:
    :caption: User Guide
    :glob:

    installation
    usage
    modules

..  toctree::
    :maxdepth: 3
    :titlesonly:
    :caption: API Reference
    :glob:

    api/modules

Indices and tables
##################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. automodule:: georges

.. |Actions Status| image:: https://github.com/ULB-Metronu/georges/actions/workflows/ci.yml/badge.svg?branch=master
   :target: https://github.com/ULB-Metronu/georges/actions

.. |Documentation Status| image:: https://github.com/ULB-Metronu/georges/actions/workflows/documentation.yml/badge.svg?branch=master
   :target: https://github.com/ULB-Metronu/georges/actions

.. |Python version| image:: _static/python_versions.svg

.. |version| image:: https://img.shields.io/badge/version-2023.1-blue

.. |Bugs| image:: https://sonarcloud.io/api/project_badges/measure?project=ULB-Metronu_georges&metric=bugs
   :target: https://sonarcloud.io/summary/overall?id=ULB-Metronu_georges

.. |Coverage| image:: https://sonarcloud.io/api/project_badges/measure?project=ULB-Metronu_georges&metric=coverage
   :target: https://sonarcloud.io/summary/overall?id=ULB-Metronu_georges

.. |Reliability| image:: https://sonarcloud.io/api/project_badges/measure?project=ULB-Metronu_georges&metric=reliability_rating
   :target: https://sonarcloud.io/summary/overall?id=ULB-Metronu_georges

.. |License| image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

.. |Code Style| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black

.. |Gitter| image:: https://badges.gitter.im/ULB-Metronu/georges.svg?
   :target: https://gitter.im/ULB-Metronu/georges?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge

