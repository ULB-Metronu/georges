*******
Manzoni
*******

The Manzoni module aims to implement fast particle beam tracking through the most encountered accelerator
and beamline elements (magnets, scatterers and degraders, collimators and cavities).

Structure
=========

The module is structured as depicted in figure below and described in the next sections:


.. figure:: ../_static/georges_structure.png
    :align: center
    :width: 100%

* :ref:`elements`
* :ref:`apertures`
* :ref:`beam`
* :ref:`input`
* :ref:`core`
* :ref:`integrators`
* :ref:`kernels`
* :ref:`maps`
* :ref:`observers`

..  toctree::
    :hidden:
    :maxdepth: 1

    manzoni/structure


Examples
========

..  toctree::
    :maxdepth: 1

    manzoni/example_tracking
    manzoni/example_energy_degradation

Validation
==========

..  toctree::
    :maxdepth: 1

    manzoni/madx_validation