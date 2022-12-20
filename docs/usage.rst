Using and executing Georges
===========================

You can access the library by simply importing it:

.. jupyter-execute::
    :hide-output:

    import georges

This will include only the core components of Georges. When importing `georges`, it will automatically import the
`georges.Kinematics` and `georges.Distribution` modules. These modules are described in details in
the Georges-core's documentation. You can simply use these modules by typing:

.. jupyter-execute::
    :hide-output:

    _ureg = georges.ureg
    kin = georges.Kinematics(120 * _ureg.MeV)
    dist = georges.Distribution()

Depending on your needs, the different Georges' modules must be imported separately:

.. jupyter-execute::
    :hide-output:

    import georges.fermi
    import georges.manzoni
    import georges.vis
