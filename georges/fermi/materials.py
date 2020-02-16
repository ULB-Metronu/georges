"""
TODO
"""
from typing import List, Mapping
import os
import sys
import numpy as _np
import pandas as _pd
from scipy.integrate import quad as _quad
import scipy.interpolate
from .mcs import ScatteringModelType as _ScatteringModelType
from .mcs import DifferentialMoliere as _DifferentialMoliere
from georges_core import ureg as _ureg
from georges_core import Kinematics as _Kinematics
from georges_core.kinematics import ekin_to_pv as _ekin_to_pv


def __pdg_read_data(path: str) -> _pd.DataFrame:
    return _pd.read_csv(os.path.join(path, "pdg", "data.csv"),
                        delimiter=',',
                        index_col='material'
                        )


def __bdsim_read_data(path: str) -> _pd.DataFrame:
    return _pd.read_csv(os.path.join(path, "bdsim", "data.csv"),
                        delimiter=',',
                        index_col='material'
                        )


__STAR_COLUMNS: List[str] = ['K', 'eS', 'nS', 'tS', 'csda', 'prange', 'factor']


def __star_read_data(path: str, material_name: str) -> _pd.DataFrame:
    return _pd.read_csv(os.path.join(path, "pstar", f"{material_name}.txt"),
                        skiprows=7,
                        delimiter=' ',
                        names=__STAR_COLUMNS,
                        index_col=False)


def __srim_read_data(path: str, material_name: str, pdg_data) -> _pd.DataFrame:
    data = _pd.read_csv(
        os.path.join(path, "srim", f"{material_name}.txt"),
        skiprows=23,
        skipfooter=20,
        engine='python',
        delim_whitespace=True,
        names=[
            'energy',
            'energy_unit',
            'stopping_elec',
            'stopping_nuclear',
            'projected_range',
            'projected_range_unit',
            'longitudinal_straggling',
            'longitudinal_straggling_unit',
            'lateral_straggling',
            'lateral_straggling_unit',
        ]
    )
    units = {
        'A': 1e-10,
        'um': 1e-6,
        'mm': 1e-3,
        'm': 1,
        'keV': 1e-3,
        'MeV': 1,
        'GeV': 1000,
    }

    data['projected_range'] = data['projected_range'].astype(float)
    data['projected_range_meter'] = data.apply(
        lambda e: units[e['projected_range_unit']] * pdg_data.at[material_name, 'rho'] * 100 * e['projected_range'],
        axis=1
    )
    data['energy'] = data['energy'].astype(float)
    data['energy_scale'] = data.apply(lambda e: units[e['energy_unit']] * e['energy'], axis=1)
    data['stopping_elec'] = data['stopping_elec'].astype(float)
    data['stopping_nuclear'] = data['stopping_nuclear'].astype(float)

    data['stopping_elec'] *= 1000
    data['stopping_nuclear'] *= 1000
    return data


class MaterialType(type):
    @property
    def valid_data(cls):
        return cls.pdg_data is not None and cls.projected_range is not None # and cls.csda_range is not None

    @property
    def atomic_a(cls) -> float:
        if cls.pdg_data is None:
            return None
        return cls.pdg_data['A']

    @property
    def atomic_z(cls) -> float:
        if cls.pdg_data is None:
            return None
        return cls.pdg_data['Z']

    @property
    def density(cls) -> float:
        if cls.pdg_data is None:
            return None
        return cls.pdg_data['rho']

    def range(cls, kinetic_energy: _ureg.Quantity, csda: bool = False, projected: bool = True) -> _ureg.Quantity:
        """

        Args:
            kinetic_energy:
            csda:
            projected:

        Returns:

        """
        if not cls.valid_data:
            return None
        energy = kinetic_energy.m_as('MeV')
        if projected and not csda:
            return _np.exp(cls.projected_range(_np.log(energy))) / cls.density * _ureg.cm
        elif csda and not projected:
            return _np.exp(cls.csda_range(_np.log(energy))) / cls.density * _ureg.cm
        else:
            raise Exception("'projected' or 'csda' arguments are mutually exclusive and one must be defined.")

    def solve_range(cls, range: _ureg.Quantity, csda: bool = False, projected: bool = True) -> _Kinematics:
        """

        Args:
            range:
            csda:
            projected:

        Returns:

        """
        if not cls.valid_data:
            return None
        r = range.m_as('cm')
        if projected and not csda:
            return _Kinematics(
                _np.exp(cls.projected_range.solve(_np.log(r * cls.density), extrapolate=False))[0] * _ureg.MeV,
                kinetic=True
            )
        elif csda and not projected:
            return _Kinematics(
                _np.exp(cls.csda_range.solve(_np.log(r * cls.density), extrapolate=False))[0] * _ureg.MeV,
                kinetic=True
            )
        else:
            raise Exception("'projected' or 'csda' arguments are mutually exclusive and one must be defined.")

    def stopping(cls, thickness: _ureg.Quantity, kinetic_energy: _ureg.Quantity) -> _Kinematics:
        """

        Args:
            thickness:
            kinetic_energy:

        Returns:

        """
        if not cls.valid_data:
            return None
        return cls.solve_range(cls.residual_range(thickness, kinetic_energy))

    def residual_range(cls, thickness: _ureg.Quantity, kinetic_energy: _ureg.Quantity) -> _ureg.Quantity:
        """

        Args:
            thickness:
            kinetic_energy:

        Returns:

        """
        if not cls.valid_data:
            return None
        return cls.range(kinetic_energy) - thickness

    def required_thickness(cls,
                           kinetic_energy_out: _ureg.Quantity,
                           kinetic_energy_in: _ureg.Quantity,
                           csda: bool = False,
                           projected: bool = True):
        """

        Args:
            kinetic_energy_out:
            kinetic_energy_in:
            csda:
            projected:

        Returns:

        """
        if not cls.valid_data:
            return None
        return cls.range(kinetic_energy_in, csda=csda, projected=projected) \
               - cls.range(kinetic_energy_out, csda=csda, projected=projected)

    def scattering(cls,
                   kinetic_energy: _ureg.Quantity,
                   thickness: _ureg.Quantity,
                   model: _ScatteringModelType = _DifferentialMoliere) -> Mapping[str, float]:
        """
        Compute the Fermi-Eyges parameters A0, A1, A2 and B (emittance).

        Args:
            kinetic_energy:
            thickness:
            model:

        Returns:

        """
        if not cls.valid_data:
            """ TODO """
            "Set a warning message"
            return {
                'A': (0.0, 0.0, 0.0),
                'B': 0.0,
                'TWISS_ALPHA': _np.nan,
                'TWISS_BETA': _np.nan,
                'TWISS_GAMMA': _np.nan,
            }
        thickness = thickness.m_as('cm')

        def integrand(u: float,
                      initial_energy: _ureg.Quantity,
                      thickness: float,
                      material: MaterialType,
                      scattering_model: _ScatteringModelType,
                      n: int):
            return (thickness - u) ** n * scattering_model.t(
                _ekin_to_pv(cls.stopping(u * _ureg.cm, initial_energy).ekin).m_as('MeV'),
                _ekin_to_pv(initial_energy).m_as('MeV'),
                material=material
            )

        a = [
            _quad(integrand, 0, thickness, args=(kinetic_energy, thickness, cls, model, 0))[0],  # Order 0
            1e-2 * _quad(integrand, 0, thickness, args=(kinetic_energy, thickness, cls, model, 1))[0],  # Order 1
            1e-4 * _quad(integrand, 0, thickness, args=(kinetic_energy, thickness, cls, model, 2))[0],  # Order 2
        ]
        b = _np.sqrt(a[0] * a[2] - a[1] ** 2)  # Emittance in m.rad

        return {
            'A': a,
            'B': b,
            'TWISS_ALPHA': -a[1] / b,
            'TWISS_BETA': a[2] / b,
            'TWISS_GAMMA': a[0] / b,
        }

    def energy_dispersion(cls, energy):
        """
        TODO

        Args:
            energy:

        Returns:

        """
        if not cls._energy_dispersion:
            return None
        c0 = cls.bdsim_data['C0']
        c1 = cls.bdsim_data['C1']
        c2 = cls.bdsim_data['C2']
        c3 = cls.bdsim_data['C3']
        c4 = cls.bdsim_data['C4']
        return c4 * energy**4 + c3 * energy**3 + c2 * energy**2 + c1 * energy + c0

    def losses(cls, energy):
        """
        TODO

        Args:
            energy:

        Returns:

        """
        if not cls._loss_factors:
            return None
        c0 = cls.bdsim_data['C0']
        c1 = cls.bdsim_data['C1']
        c2 = cls.bdsim_data['C2']
        c3 = cls.bdsim_data['C3']
        c4 = cls.bdsim_data['C4']
        return c4 * energy ** 4 + c3 * energy ** 3 + c2 * energy ** 2 + c1 * energy + c0


class Material(metaclass=MaterialType):
    def __init_subclass__(cls,
                          csda_range_data,
                          projected_range_data,
                          pdg_data,
                          bdsim_data,
                          **kwargs):
        super().__init_subclass__(**kwargs)
        cls.csda_range = csda_range_data
        cls.projected_range = projected_range_data
        cls.pdg_data = pdg_data
        cls.bdsim_data = bdsim_data

    def __str__(self):
        return self.__class__.__name__.lower()

    def __eq__(self, y):
        return str(self) == str(y)


def __initialize_materials_database():
    db_path = os.path.dirname(__file__)

    # Read PDG data
    pdg_data = __pdg_read_data(db_path)

    # Read P-Star data
    pstar_data_files = [os.path.splitext(f)[0] for f in os.listdir(os.path.join(db_path, 'pstar')) if
                        os.path.isfile(os.path.join(db_path, 'pstar', f)) and os.path.splitext(f)[1] == '.txt']
    pstar = {
        m: __star_read_data(db_path, m) for m in pstar_data_files
    }

    # Read SRIM data
    srim_data_files = [os.path.splitext(f)[0] for f in os.listdir(os.path.join(db_path, 'srim')) if
                       os.path.isfile(os.path.join(db_path, 'srim', f)) and os.path.splitext(f)[1] == '.txt']
    srim = {
        m: __srim_read_data(db_path, m, pdg_data) for m in srim_data_files
    }

    # Read BDSIM data
    bdsim_data = __bdsim_read_data(db_path)

    # Create the interpolated functions for the projected ranges
    projected_ranges_pstar = {
        k: scipy.interpolate.CubicSpline(_np.log(v['K']), _np.log(v['prange'])) for k, v in
        pstar.items()
    }
    projected_ranges_srim = {
        k: scipy.interpolate.CubicSpline(_np.log(v['energy_scale']), _np.log(v['projected_range_meter'])) for k, v in
        srim.items()
    }
    projected_ranges = {**projected_ranges_pstar, **projected_ranges_srim}

    # Create the interpolated functions for the CSDA ranges
    csda_ranges_pstar = {
        k: scipy.interpolate.CubicSpline(_np.log(v['K']), _np.log(v['csda'])) for k, v in
        pstar.items()
    }
    csda_ranges = {**csda_ranges_pstar}

    # Dynamically create the material classes
    for m in set(('Vacuum',) + tuple(projected_ranges.keys()) + tuple(csda_ranges.keys())):
        setattr(sys.modules[__name__], m.title(), MaterialType(m.title(), (Material,), {},
                                                               projected_range_data=projected_ranges.get(m, None),
                                                               csda_range_data=csda_ranges.get(m, None),
                                                               pdg_data=pdg_data.loc[
                                                                   m] if m in pdg_data.index else None,
                                                               bdsim_data=bdsim_data.loc[
                                                                   m] if m in bdsim_data.index else None
                                                               ))
