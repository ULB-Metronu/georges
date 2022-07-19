"""
Testing procedure of the Fermi module.
Data are coming from the paper: On the scattering Power of radiotherapy protons from B.Gottschalk
"""
import pytest
import numpy as _np
import georges
import georges.fermi
import georges.fermi.materials as gmat
import georges.fermi.mcs as gmod
from georges import ureg as _ureg


@pytest.mark.parametrize("material, model, ratio, epos_expected, angle_hanson, angle_model", [
    # FermiRossi
    (gmat.Beryllium, gmod.FermiRossi, 0.001, 158.51, 0.565, 62.86),
    (gmat.Beryllium, gmod.FermiRossi, 0.010, 157.69, 2.020, 44.55),
    (gmat.Beryllium, gmod.FermiRossi, 0.100, 149.32, 7.207, 31.43),
    (gmat.Beryllium, gmod.FermiRossi, 0.200, 139.62, 10.784, 28.12),
    (gmat.Beryllium, gmod.FermiRossi, 0.500, 107.00, 19.812, 23.93),
    (gmat.Beryllium, gmod.FermiRossi, 0.900, 43.71, 37.994, 21.34),
    (gmat.Beryllium, gmod.FermiRossi, 0.970, 22.49, 48.051, 21.33),
    # ICRUProtons
    (gmat.Beryllium, gmod.ICRUProtons, 0.001, 158.51, 0.565, 36.65),
    (gmat.Beryllium, gmod.ICRUProtons, 0.010, 157.69, 2.020, 21.29),
    (gmat.Beryllium, gmod.ICRUProtons, 0.100, 149.32, 7.207, 10.28),
    (gmat.Beryllium, gmod.ICRUProtons, 0.200, 139.62, 10.784, 7.50),
    (gmat.Beryllium, gmod.ICRUProtons, 0.500, 107.00, 19.812, 3.98),
    (gmat.Beryllium, gmod.ICRUProtons, 0.900, 43.71, 37.994, 1.82),
    (gmat.Beryllium, gmod.ICRUProtons, 0.970, 22.49, 48.051, 1.80),
    # DiffHighland
    (gmat.Beryllium, gmod.DifferentialHighland, 0.001, 158.51, 0.565, -6.51),
    (gmat.Beryllium, gmod.DifferentialHighland, 0.010, 157.69, 2.020, -1.89),
    (gmat.Beryllium, gmod.DifferentialHighland, 0.100, 149.32, 7.207, 3.10),
    (gmat.Beryllium, gmod.DifferentialHighland, 0.200, 139.62, 10.784, 4.72),
    (gmat.Beryllium, gmod.DifferentialHighland, 0.500, 107.00, 19.812, 7.18),
    (gmat.Beryllium, gmod.DifferentialHighland, 0.900, 43.71, 37.994, 9.81),
    (gmat.Beryllium, gmod.DifferentialHighland, 0.970, 22.49, 48.051, 11.0),
    # DiffMoliere
    (gmat.Beryllium, gmod.DifferentialMoliere, 0.001, 158.51, 0.565, -3.03),
    (gmat.Beryllium, gmod.DifferentialMoliere, 0.010, 157.69, 2.020, -0.22),
    (gmat.Beryllium, gmod.DifferentialMoliere, 0.100, 149.32, 7.207, 1.63),
    (gmat.Beryllium, gmod.DifferentialMoliere, 0.200, 139.62, 10.784, 2.02),
    (gmat.Beryllium, gmod.DifferentialMoliere, 0.500, 107.00, 19.812, 2.21),
    (gmat.Beryllium, gmod.DifferentialMoliere, 0.900, 43.71, 37.994, 1.14),
    (gmat.Beryllium, gmod.DifferentialMoliere, 0.970, 22.49, 48.051, 0.14),
    # FermiRossi
    (gmat.Aluminum, gmod.FermiRossi, 0.001, 158.51, 1.004, 54.97),
    (gmat.Aluminum, gmod.FermiRossi, 0.010, 157.68, 3.568, 34.83),
    (gmat.Aluminum, gmod.FermiRossi, 0.100, 149.22, 13.214, 21.11),
    (gmat.Aluminum, gmod.FermiRossi, 0.200, 139.42, 19.836, 17.72),
    (gmat.Aluminum, gmod.FermiRossi, 0.500, 106.52, 36.631, 13.45),
    (gmat.Aluminum, gmod.FermiRossi, 0.900, 42.99, 70.843, 10.78),
    (gmat.Aluminum, gmod.FermiRossi, 0.970, 21.83, 90.068, 10.81),
    # ICRUProtons
    (gmat.Aluminum, gmod.ICRUProtons, 0.001, 158.51, 1.004, 41.63),
    (gmat.Aluminum, gmod.ICRUProtons, 0.010, 157.68, 3.568, 23.22),
    (gmat.Aluminum, gmod.ICRUProtons, 0.100, 149.22, 13.214, 10.69),
    (gmat.Aluminum, gmod.ICRUProtons, 0.200, 139.42, 19.836, 7.59),
    (gmat.Aluminum, gmod.ICRUProtons, 0.500, 106.52, 36.631, 3.68),
    (gmat.Aluminum, gmod.ICRUProtons, 0.900, 42.99, 70.843, 1.25),
    (gmat.Aluminum, gmod.ICRUProtons, 0.970, 21.83, 90.068, 1.27),
    # DiffHighland
    (gmat.Aluminum, gmod.DifferentialHighland, 0.001, 158.51, 1.004, -3.67),
    (gmat.Aluminum, gmod.DifferentialHighland, 0.010, 157.68, 3.568, -2.08),
    (gmat.Aluminum, gmod.DifferentialHighland, 0.100, 149.22, 13.214, 0.77),
    (gmat.Aluminum, gmod.DifferentialHighland, 0.200, 139.42, 19.836, 1.83),
    (gmat.Aluminum, gmod.DifferentialHighland, 0.500, 106.52, 36.631, 3.53),
    (gmat.Aluminum, gmod.DifferentialHighland, 0.900, 42.99, 70.843, 5.56),
    (gmat.Aluminum, gmod.DifferentialHighland, 0.970, 21.83, 90.068, 6.70),
    # DiffMoliere
    (gmat.Aluminum, gmod.DifferentialMoliere, 0.001, 158.51, 1.004, 0.56),
    (gmat.Aluminum, gmod.DifferentialMoliere, 0.010, 157.68, 3.568, 1.43),
    (gmat.Aluminum, gmod.DifferentialMoliere, 0.100, 149.22, 13.214, 2.06),
    (gmat.Aluminum, gmod.DifferentialMoliere, 0.200, 139.42, 19.836, 2.14),
    (gmat.Aluminum, gmod.DifferentialMoliere, 0.500, 106.52, 36.631, 1.94),
    (gmat.Aluminum, gmod.DifferentialMoliere, 0.900, 42.99, 70.843, 0.56),
    (gmat.Aluminum, gmod.DifferentialMoliere, 0.970, 21.83, 90.068, -0.45),
    # FermiRossi
    (gmat.Copper, gmod.FermiRossi, 0.001, 158.51, 1.545, 49.12),
    (gmat.Copper, gmod.FermiRossi, 0.010, 157.67, 5.664, 28.90),
    (gmat.Copper, gmod.FermiRossi, 0.100, 149.14, 20.535, 15.40),
    (gmat.Copper, gmod.FermiRossi, 0.200, 139.25, 30.855, 12.10),
    (gmat.Copper, gmod.FermiRossi, 0.500, 106.08, 57.059, 8.02),
    (gmat.Copper, gmod.FermiRossi, 0.900, 42.22, 110.652, 5.85),
    (gmat.Copper, gmod.FermiRossi, 0.970, 21.14, 140.987, 6.32),
    # ICRUProtons
    (gmat.Copper, gmod.ICRUProtons, 0.001, 158.51, 1.545, 39.85),
    (gmat.Copper, gmod.ICRUProtons, 0.010, 157.67, 5.664, 20.88),
    (gmat.Copper, gmod.ICRUProtons, 0.100, 149.14, 20.535, 8.22),
    (gmat.Copper, gmod.ICRUProtons, 0.200, 139.25, 30.855, 5.13),
    (gmat.Copper, gmod.ICRUProtons, 0.500, 106.08, 57.059, 1.31),
    (gmat.Copper, gmod.ICRUProtons, 0.900, 42.22, 110.652, -0.73),
    (gmat.Copper, gmod.ICRUProtons, 0.970, 21.14, 140.987, -0.29),
    # DiffHighland
    (gmat.Copper, gmod.DifferentialHighland, 0.001, 158.51, 1.545, -1.99),
    (gmat.Copper, gmod.DifferentialHighland, 0.010, 157.67, 5.664, -1.79),
    (gmat.Copper, gmod.DifferentialHighland, 0.100, 149.14, 20.535, 0.13),
    (gmat.Copper, gmod.DifferentialHighland, 0.200, 139.25, 30.855, 0.97),
    (gmat.Copper, gmod.DifferentialHighland, 0.500, 106.08, 57.059, 2.43),
    (gmat.Copper, gmod.DifferentialHighland, 0.900, 42.22, 110.652, 4.66),
    (gmat.Copper, gmod.DifferentialHighland, 0.970, 21.14, 140.987, 6.20),
    # DiffMoliere
    (gmat.Copper, gmod.DifferentialMoliere, 0.001, 158.51, 1.545, -0.60),
    (gmat.Copper, gmod.DifferentialMoliere, 0.010, 157.67, 5.664, -0.44),
    (gmat.Copper, gmod.DifferentialMoliere, 0.100, 149.14, 20.535, -0.18),
    (gmat.Copper, gmod.DifferentialMoliere, 0.200, 139.25, 30.855, -0.15),
    (gmat.Copper, gmod.DifferentialMoliere, 0.500, 106.08, 57.059, -0.36),
    (gmat.Copper, gmod.DifferentialMoliere, 0.900, 42.22, 110.652, -1.41),
    (gmat.Copper, gmod.DifferentialMoliere, 0.970, 21.14, 140.987, -2.07),
    # FermiRossi
    (gmat.Lead, gmod.FermiRossi, 0.001, 158.51, 2.638, 45.33),
    (gmat.Lead, gmod.FermiRossi, 0.010, 157.65, 9.878, 23.07),
    (gmat.Lead, gmod.FermiRossi, 0.100, 148.96, 36.222, 8.99),
    (gmat.Lead, gmod.FermiRossi, 0.200, 138.91, 54.568, 5.66),
    (gmat.Lead, gmod.FermiRossi, 0.500, 105.27, 101.274, 1.71),
    (gmat.Lead, gmod.FermiRossi, 0.900, 40.97, 197.440, 0.15),
    (gmat.Lead, gmod.FermiRossi, 0.970, 19.95, 253.086, 1.15),
    # ICRUProtons
    (gmat.Lead, gmod.ICRUProtons, 0.001, 158.51, 2.638, 42.50),
    (gmat.Lead, gmod.ICRUProtons, 0.010, 157.65, 9.878, 20.68),
    (gmat.Lead, gmod.ICRUProtons, 0.100, 148.96, 36.222, 6.87),
    (gmat.Lead, gmod.ICRUProtons, 0.200, 138.91, 54.568, 3.60),
    (gmat.Lead, gmod.ICRUProtons, 0.500, 105.27, 101.274, -0.27),
    (gmat.Lead, gmod.ICRUProtons, 0.900, 40.97, 197.440, -1.80),
    (gmat.Lead, gmod.ICRUProtons, 0.970, 19.95, 253.086, -0.82),
    # DiffHighland
    (gmat.Lead, gmod.DifferentialHighland, 0.001, 158.51, 2.638, -2.24),
    (gmat.Lead, gmod.DifferentialHighland, 0.010, 157.65, 9.878, -0.53),
    (gmat.Lead, gmod.DifferentialHighland, 0.100, 148.96, 36.222, -0.38),
    (gmat.Lead, gmod.DifferentialHighland, 0.200, 138.91, 54.568, 0.06),
    (gmat.Lead, gmod.DifferentialHighland, 0.500, 105.27, 101.274, 1.17),
    (gmat.Lead, gmod.DifferentialHighland, 0.900, 40.97, 197.440, 3.71),
    (gmat.Lead, gmod.DifferentialHighland, 0.970, 19.95, 253.086, 5.78),
    # DiffMoliere
    (gmat.Lead, gmod.DifferentialMoliere, 0.001, 158.51, 2.638, 1.39),
    (gmat.Lead, gmod.DifferentialMoliere, 0.010, 157.65, 9.878, -0.51),
    (gmat.Lead, gmod.DifferentialMoliere, 0.100, 148.96, 36.222, -1.35),
    (gmat.Lead, gmod.DifferentialMoliere, 0.200, 138.91, 54.568, -1.53),
    (gmat.Lead, gmod.DifferentialMoliere, 0.500, 105.27, 101.274, -1.86),
    (gmat.Lead, gmod.DifferentialMoliere, 0.900, 40.97, 197.440, -2.51),
    (gmat.Lead, gmod.DifferentialMoliere, 0.970, 19.95, 253.086, -2.74),
])
def test_angle(material, model, ratio, epos_expected, angle_hanson, angle_model):
    kin = georges.Kinematics(158.6 * _ureg.MeV)
    thickness = ratio * material.range(kinetic_energy=kin.ekin)
    epos = material.stopping(thickness=thickness, kinetic_energy=kin.ekin).ekin.m_as('MeV')

    angle = material.scattering(kinetic_energy=kin.ekin, thickness=thickness, model=model)
    a2 = 1000 * (_np.sqrt(angle['A'][0]))

    _np.testing.assert_allclose(epos, epos_expected, rtol=1e-2)
    _np.testing.assert_allclose(a2, angle_hanson * ((angle_model / 100) + 1), rtol=5e-2)
