"""
TODO
"""
from typing import Optional, List, Mapping
import os
import sys
import yaml
import numpy as _np
import pandas as _pd
import scipy.constants
from scipy.integrate import quad as _quad
import scipy.interpolate
from .mcs import ScatteringModelType as _ScatteringModelType
from .mcs import DifferentialMoliere as _DifferentialMoliere
from georges_core import ureg as _ureg
from georges_core import Kinematics as _Kinematics
from georges_core.kinematics import ekin_to_pv as _ekin_to_pv


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


def __srim_read_data(path: str, material_name: str, density: float) -> _pd.DataFrame:
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
        lambda e: units[e['projected_range_unit']] * density * 100 * e['projected_range'],
        axis=1
    )
    data['energy'] = data['energy'].astype(float)
    data['energy_scale'] = data.apply(lambda e: units[e['energy_unit']] * e['energy'], axis=1)
    data['stopping_elec'] = data['stopping_elec'].astype(float)
    data['stopping_nuclear'] = data['stopping_nuclear'].astype(float)

    data['stopping_elec'] *= 1000
    data['stopping_nuclear'] *= 1000
    return data


class RangeDefinitionType(type):
    pass


class CSDARange(metaclass=RangeDefinitionType):
    pass


class ProjectedRange(metaclass=RangeDefinitionType):
    pass


class ElementType(type):
    @property
    def atomic_a(cls) -> float:
        if cls.material_data is None:
            return None
        return cls.material_data['A']

    @property
    def atomic_z(cls) -> int:
        if cls.material_data is None:
            return None
        return cls.material_data['Z']

    @property
    def atomic_radiation_length(cls) -> float:
        """

        Args:
            material:

        Returns:

        """
        a: float = cls.atomic_a
        z: int = cls.atomic_z
        # Probably wrong, should include rho somewhere
        return 716.4 * a * (1 / (z * (z + 1) * (_np.log(287 / _np.sqrt(z)))))

    @property
    def atomic_scattering_length(cls) -> float:
        """
        See "Techniques of Proton Radiotherapy:Transport Theory", B. Gottschalk, 2012.

        Args:
            material:

        Returns:

        """
        alpha: float = scipy.constants.alpha  # Fine structure constant
        avogadro: float = scipy.constants.Avogadro
        re: float = scipy.constants.physical_constants['classical electron radius'][0] * 100  # in cm!
        a: float = cls.atomic_a
        z: int = cls.atomic_z
        return 1 / (alpha * avogadro * re ** 2 * z ** 2 * (2 * _np.log(33219 * (a * z) ** (-1 / 3)) - 1) / a)


class Element(metaclass=ElementType):
    def __init_subclass__(cls,
                          material_data,
                          **kwargs):
        super().__init_subclass__(**kwargs)
        cls.material_data = material_data

    def __str__(self):
        return self.__class__.__name__.lower()

    def __eq__(self, y):
        return str(self) == str(y)


class CompoundType(type):
    @property
    def valid_data(cls):
        return cls.material_data is not None and cls.projected_range is not None and (cls.projected_range is not None or cls.csda_range is not None)

    @property
    def density(cls) -> float:
        if cls.material_data is None:
            return None
        return cls.material_data['rho']

    @density.setter
    def density(cls, value) -> float:
        cls.material_data['rho'] = value

    @property
    def radiation_length(cls) -> float:
        length = 0
        for element, fraction in cls.material_data['fractions'].items():
            element = getattr(sys.modules[__name__], 'elements')[element.title()]
            length += 1.0 / fraction * element.atomic_radiation_length
        return 1/length

    @property
    def scattering_length(cls) -> float:
        length = 0
        for element, fraction in cls.material_data['fractions'].items():
            element = getattr(sys.modules[__name__], 'elements')[element.title()]
            length += fraction / element.atomic_scattering_length
        return 1/(cls.density * length)

    def range(cls,
              kinetic_energy: _ureg.Quantity,
              range_definition: RangeDefinitionType = ProjectedRange
              ) -> Optional[_ureg.Quantity]:
        """

        Args:
            kinetic_energy:
            range_definition:

        Returns:

        """
        if not cls.valid_data:
            return None
        energy = kinetic_energy.m_as('MeV')
        if range_definition == ProjectedRange:
            return _np.exp(cls.projected_range(_np.log(energy))) / cls.density * _ureg.cm
        elif range_definition == CSDARange:
            return _np.exp(cls.csda_range(_np.log(energy))) / cls.density * _ureg.cm
        else:
            raise Exception("'projected' or 'csda' arguments are mutually exclusive and one must be defined.")

    def solve_range(cls,
                    range: _ureg.Quantity,
                    range_definition: RangeDefinitionType = ProjectedRange
                    ) -> Optional[_Kinematics]:
        """

        Args:
            range:
            range_definition:

        Returns:

        """
        if not cls.valid_data:
            return None
        r = range.m_as('cm')
        if range_definition == ProjectedRange:
            return _Kinematics(
                _np.exp(cls.projected_range.solve(_np.log(r * cls.density), extrapolate=False))[0] * _ureg.MeV,
                kinetic=True
            )
        elif range_definition == CSDARange:
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
            return _Kinematics(kinetic_energy)
        return cls.solve_range(cls.residual_range(thickness, kinetic_energy))

    def residual_range(cls, thickness: _ureg.Quantity, kinetic_energy: _ureg.Quantity) -> Optional[_ureg.Quantity]:
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
                           range_definition: RangeDefinitionType = ProjectedRange):
        """

        Args:
            kinetic_energy_out:
            kinetic_energy_in:
            range_definition:

        Returns:

        """
        if not cls.valid_data:
            return None
        return cls.range(kinetic_energy_in, range_definition=range_definition) \
               - cls.range(kinetic_energy_out, range_definition=range_definition)

    def scattering(cls,
                   kinetic_energy: _ureg.Quantity,
                   thickness: _ureg.Quantity,
                   model: _ScatteringModelType = _DifferentialMoliere,
                   compute_a0: bool = True,
                   compute_a1: bool = True,
                   compute_a2: bool = True) -> Mapping[str, float]:
        """
        Compute the Fermi-Eyges parameters A0, A1, A2 and B (emittance).

        Args:
            kinetic_energy:
            thickness:
            model:
            compute_a0:
            compute_a1:
            compute_a2:

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
                      material: ElementType,
                      scattering_model: _ScatteringModelType,
                      n: int):
            return (thickness - u) ** n * scattering_model.t(
                _ekin_to_pv(cls.stopping(u * _ureg.cm, initial_energy).ekin).m_as('MeV'),
                _ekin_to_pv(initial_energy).m_as('MeV'),
                material=material
            )
        a = [0.0, 0.0, 0.0]
        if compute_a0:
            a[0] = _quad(integrand, 0, thickness, args=(kinetic_energy, thickness, cls, model, 0))[0]  # Order 0
        if compute_a1:
            a[1] = 1e-2 * _quad(integrand, 0, thickness, args=(kinetic_energy, thickness, cls, model, 1))[0]  # Order 1
        if compute_a2:
            a[2] = 1e-4 * _quad(integrand, 0, thickness, args=(kinetic_energy, thickness, cls, model, 2))[0]  # Order 2

        if compute_a0 and compute_a1 and compute_a2:
            b = _np.sqrt(a[0] * a[2] - a[1] ** 2)  # Emittance in m.rad
        else:
            b = 0.0

        return {
            'A': a,
            'B': b,
            'TWISS_ALPHA': -a[1] / b if b != 0.0 else 0.0,
            'TWISS_BETA': a[2] / b if b != 0.0 else 0.0,
            'TWISS_GAMMA': a[0] / b if b != 0.0 else 0.0,
        }

    def energy_dispersion(cls, energy: _ureg.Quantity) -> float:
        """
        TODO

        Args:
            energy:

        Returns:

        """
        if cls.bdsim_data is None:
            return 0.0

        c0 = cls.bdsim_data['energy_dispersion_c0']
        c1 = cls.bdsim_data['energy_dispersion_c1']
        c2 = cls.bdsim_data['energy_dispersion_c2']
        c3 = cls.bdsim_data['energy_dispersion_c3']
        c4 = cls.bdsim_data['energy_dispersion_c4']
        c5 = cls.bdsim_data['energy_dispersion_c5']
        energy = energy.m_as("MeV")
        return c5 * energy**5 + c4 * energy**4 + c3 * energy**3 + c2 * energy**2 + c1 * energy + c0

    def losses(cls, energy: _ureg.Quantity) -> float:
        """
        TODO

        Args:
            energy:

        Returns:

        """
        if cls.bdsim_data is None:
            return 1.0

        c0 = cls.bdsim_data['loss_c0']
        c1 = cls.bdsim_data['loss_c1']
        c2 = cls.bdsim_data['loss_c2']
        c3 = cls.bdsim_data['loss_c3']
        energy = energy.m_as("MeV")
        return c3 * energy**3 + c2 * energy**2 + c1 * energy + c0


class Compound(metaclass=CompoundType):
    def __init_subclass__(cls,
                          csda_range_data,
                          projected_range_data,
                          material_data,
                          bdsim_data,
                          **kwargs):
        super().__init_subclass__(**kwargs)
        cls.csda_range = csda_range_data
        cls.projected_range = projected_range_data
        cls.material_data = material_data
        cls.bdsim_data = bdsim_data

    def __str__(self):
        return self.__class__.__name__.lower()

    def __eq__(self, y):
        return str(self) == str(y)


def __initialize_materials_database():
    db_path = os.path.dirname(__file__)

    # Read all material definitions
    materials_definitions = yaml.safe_load(open(os.path.join(db_path, "materials.yaml"), 'r'))

    # Read P-Star data
    pstar_data_files = [os.path.splitext(f)[0] for f in os.listdir(os.path.join(db_path, 'pstar')) if
                        os.path.isfile(os.path.join(db_path, 'pstar', f)) and os.path.splitext(f)[1] == '.txt']
    pstar = {
        m: __star_read_data(db_path, m) for m in pstar_data_files if materials_definitions['compounds'].get(m)
    }

    # Read SRIM data
    srim_data_files = [os.path.splitext(f)[0] for f in os.listdir(os.path.join(db_path, 'srim')) if
                       os.path.isfile(os.path.join(db_path, 'srim', f)) and os.path.splitext(f)[1] == '.txt']
    srim = {
        m: __srim_read_data(db_path,
                            m,
                            materials_definitions['compounds'].get(m)['rho']
                            ) for m in srim_data_files if materials_definitions['compounds'].get(m)
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
    elements = {}
    for k, v in materials_definitions['elements'].items():
        elements[k.title()] = ElementType(k.title(), (Element,), {},
                                                              material_data=v,
                                                              )
        setattr(sys.modules[__name__], 'elements', elements)
    for k, v in materials_definitions['compounds'].items():
        setattr(sys.modules[__name__], k.title(), CompoundType(k.title(), (Compound,), {},
                                                                      projected_range_data=projected_ranges.get(k, None),
                                                                      csda_range_data=csda_ranges.get(k, None),
                                                                      material_data=v,
                                                                      bdsim_data=bdsim_data.loc[
                                                                           k] if k in bdsim_data.index else None
                                                                      ))
