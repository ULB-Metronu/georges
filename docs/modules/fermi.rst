:tocdepth: 2

*****
Fermi
*****

Fermi-Eyges
###########

The `Fermi-Eyges` module is a reimplementation of the Fermi-Eyges transport framework, largely based on the work and publications of Bernard Gottschalk.
See :cite:p:`Hernalsteens2022` for the complete validation of this implementation.

The module is composed of:
 - *Materials database* (`materials_db.py`): reads, loads and provides an interface to data (stopping power, radiation lengths, etc.) for a large set of materials commonly used in protontherapy. The database reads stopping power and range tables from the p-star database format or from SRIM software.
 - *Stopping power* (`stopping.py`): computes the various range-energy direct and inverse relationships from measured range stopping power data tables;
 - *MCS* (multiple Coulomb scattering) (`mcs.py`): implements various scattering angle models (Gottschalk's `DifferentialMoliere` being the default);
 - *Fermi-Eyges* (`fermi_eyges.py`): computes the transport integrales (A0, A1, A2 and B) for a given material, thickness and incident energy. It also returns the residual energy to allow easy chaining through multiple slabs;
 - *Propagation* (`propagation.py`): propagation of a beam through a mutli-material beamline following the Fermi-Eyges transport theory.

Plotting support is also provided in the `georges/vis` module for the visualization of scattering beamlines.

Material database
#################

.. list-table:: Database of materials available in the Fermi Python module. The origin of the range table for each material is indicated. The density values, taken from the PDG are referenced.
    :widths: 50 25 50
    :header-rows: 1

    *   - Material
        - Range Table
        - Density (g/cm\ :sup:`3`\ )
    *   - H\ :sub:`2`\ (gaseous)
        - PSTAR
        - 0.0083748
    *   - Be
        - PSTAR
        - 1.848
    *   - B\ :sub:`4`\ C
        - SRIM
        - 2.52
    *   - Polyethylene (PE)
        - PSTAR
        - 0.94
    *   - Polystyrene (PS)
        - PSTAR
        - 1.06
    *   - C (graphite)
        - PSTAR
        - 1.7
    *   - C (diamond)
        - SRIM
        - 3.52
    *   - Polycarbonates (lexan)
        - PSTAR
        - 1.2
    *   - PET (mylar)
        - SRIM
        - 1.4
    *   - Air (gaseous)
        - PSTAR
        - 1.205e-3
    *   - Water
        - PSTAR
        - 1
    *   - O\ :sub:`2`\
        - PSTAR
        - 0.00133151
    *   - Al
        - PSTAR
        - 2.6989
    *   - Ti
        - PSTAR
        - 4.54
    *   - Co
        - PSTAR
        - 8.96
    *   - Sn
        - PSTAR
        - 7.31
    *   - Ta
        - SRIM
        - 16.66
    *   - Au
        - PSTAR
        - 19.32
    *   - Pb
        - PSTAR
        - 11.34


Usage
#####

Import the necessary submodules:

.. jupyter-execute::
   :hide-output:

    from georges.fermi import materials as gfmaterial
    from georges import ureg as _ureg

Define some materials
"""""""""""""""""""""

.. jupyter-execute::
   :hide-output:

    mat_water = gfmaterial.Water
    mat_beryllium = gfmaterial.Beryllium

Get the required thickness to degrade from 230 MeV to 70 MeV
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    thickness_be = mat_beryllium.required_thickness(70 * _ureg.MeV, 230 * _ureg.MeV)
    thickness_water = mat_water.required_thickness(70 * _ureg.MeV, 230 * _ureg.MeV)
    print(thickness_be, thickness_water)

Compute the residual kinematics after certain thickness of material
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    kpos_be = mat_beryllium.stopping(thickness=10*_ureg.cm, kinetic_energy=230*_ureg.MeV)
    kpos_water = mat_water.stopping(thickness=20*_ureg.cm, kinetic_energy=230*_ureg.MeV)
    print(kpos_be, kpos_water)

Obtain the range at a given energy
""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    print(mat_water.range(230*_ureg.MeV))
    print(mat_water.range(70*_ureg.MeV))

Obtain the kinematics required to reach a given range
"""""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    degraded_kinematics = mat_water.solve_range(15 * _ureg.cm)
    print(degraded_kinematics)


Fermi-Eyges transport theory
""""""""""""""""""""""""""""

.. jupyter-execute::

    print(mat_water.scattering(230 * _ureg.MeV, 10*_ureg.cm))
    print(mat_beryllium.scattering(100 * _ureg.MeV, 5*_ureg.cm))


Tracking in a line
******************

Let's define a line with several degraders and scatterers and we compute the parameters
`A_0`, `A_1` and `A_2` along the line.

Download: :jupyter-download:nb:`click to download <fermi>`

.. jupyter-execute::

    %matplotlib inline
    import georges
    from georges.fermi import materials
    from georges import ureg as _ureg
    from georges.manzoni.elements import Degrader

    sequence = georges.PlacementSequence(name="LINE")
    d1 = georges.Element.Degrader(NAME="D1",
                                  L=10*_ureg.cm,
                                  MATERIAL=materials.Beryllium,
                                  WITH_LOSSES=True)
    d2 = georges.Element.Scatterer(NAME="D2",
                                   L=0.1*_ureg.cm,
                                   MATERIAL=materials.Graphite)

    d3 = georges.Element.Degrader(NAME="D3",
                                  L=5*_ureg.cm,
                                  MATERIAL=materials.Aluminum,
                                  WITH_LOSSES=True)

    sequence.place(d1, at_entry=0*_ureg.m)
    sequence.place(d2, at_entry=0.5*_ureg.m)
    sequence.place(d3, at_entry=0.7*_ureg.m)

    pbs = georges.fermi.propagate(
                            sequence=sequence,
                            energy=300 *_ureg.MeV,
                            beam={
                                'A0': 0,
                                'A1': 0,
                                'A2': 0,
                            })

    s = []
    a0 = []
    a1 = []
    a2 = []
    for name, k in pbs.iterrows():
        s.append(k['AT_ENTRY'].m_as('m'))
        s.append(k['AT_EXIT'].m_as('m'))
        a0.append(k['A0_IN'])
        a0.append(k['A0_OUT'])
        a1.append(k['A1_IN'])
        a1.append(k['A1_OUT'])
        a2.append(k['A2_IN'])
        a2.append(k['A2_OUT'])

    artist = georges.vis.ManzoniMatplotlibArtist()
    artist.plot_cartouche(beamline=sequence.df)
    artist.plot(s,a0,label='a0')
    artist.plot(s,a1,label='a1')
    artist.plot(s,a2,label='a2')
    artist.ax.legend()

We can also plot the energy degradation along the line:

.. jupyter-execute::

    s = []
    edep = []
    for name, k in pbs.iterrows():
        s.append(k['AT_ENTRY'].m_as('m'))
        s.append(k['AT_EXIT'].m_as('m'))
        edep.append(k['ENERGY_IN'])
        edep.append(k['ENERGY_OUT'])

    artist = georges.vis.ManzoniMatplotlibArtist()
    artist.plot_cartouche(beamline=sequence.df)
    artist.plot(s,edep)


Python script
#############

If you would like to compute the coefficients for another material,
you must adapt the file `degrader_properties.gmad` and run the script `bdsim-input.gmad`:

::

    bdsim --file=bdsim-input.gmad --outfile=output.root --ngenerate=nparticles --batch

The program that computes the coefficients for losses and momentum
deviation is `compute_quantiles.py` and it can be excecuted by:

::

    python compute_coefficients.py path_results nparticles

Where `path_to_results` is the path to the `bdsim` output files and `nparticles` is
the number of primary particles used in the simulation.

API
###

.. automodapi:: georges.fermi
    :no-heading: