Manzoni
#######
Manzoni is a Python module that implements a matrix formalism to track proton beams through beamline magnets with
midplane symmetric magnetic fields. For each magnetic element defined in table 1, Manzoni computes the particle output
coordinates by simply applying the transfer matrix of the magnet (defined by its parameters) to the input coordinates.
This matrix formalism is made faster through the JIT compilation and is sufficiently accurate when the beam is close
to the reference (central) trajectory. Nevertheless, in general, the output coordinate does not depend linearly on
the coordinate of the incoming particle. This is typically the case in modern proton therapy techniques where the
beam is used to scan a patient’s tumor over several centimeters. Higher-order terms must then be added to the matrix
formalism of Manzoni for more precise modeling. For now, the MAD-X, the MAD-8 and the TRANSPORT integrators are
implemented.

.. csv-table:: Elements defined in Manzoni and integrators
   :header: "Element", "MAD-X", "MAD-8", "TRANSPORT"
   :widths: 40, 40, 40, 40

    "Drift", `Yes`, `Yes`, `Yes`
    "Combined Dipole", `Yes`, `Yes`, `Yes`
    "Quadrupole", `Yes`,`Yes`,`Yes`
    "Sextupole", `Yes`, `No`, `Yes`

The user should be careful to the definition of the beam. In the case of MADX,
the coordinates are *[x, px, y, py, dpp, pt]* and for MAD-8 and TRANSPORT, these are *[x, px, y, py, l, pt]*