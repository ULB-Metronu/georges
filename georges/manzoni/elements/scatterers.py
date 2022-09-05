"""
TODO
"""
from typing import Tuple, List
import numpy as _np

try:
    import numpy.random_intel as nprandom
except ModuleNotFoundError:
    import numpy.random as nprandom
from ... import ureg as _ureg
from .elements import ManzoniElement as _ManzoniElement
from ...fermi import materials
from ... import Kinematics
from ..integrators import MadXIntegrator, MadXParaxialDriftIntegrator
from ..beam import Beam as _Beam


class MaterialElement(_ManzoniElement):
    INTEGRATOR = None

    @property
    def material(self):
        if isinstance(self.MATERIAL, str):
            return getattr(materials, self.MATERIAL.capitalize())
        return self.MATERIAL

    @property
    def degraded_energy(self):
        if self.L.magnitude == 0.0:
            return self.KINETIC_ENERGY
        else:
            return (self.material.stopping(self.L, self.KINETIC_ENERGY)).ekin

    @property
    def cache(self) -> list:
        if not self.frozen:
            self._cache = self.parameters
        return self._cache

    @property
    def beta(self):
        return Kinematics(self.KINETIC_ENERGY, kinetic=True).beta

    @staticmethod
    def pseudo_aperture_check(beam, beta):
        px = beam[:, 1]
        py = beam[:, 3]
        pt = beam[:, 5]
        condition = 1 + 2 * pt / beta + pt ** 2 - px ** 2 - py ** 2 >= 0
        checked_beam = _np.compress(condition, beam, axis=0)

        return checked_beam


class Scatterer(MaterialElement):
    PARAMETERS = {
        'MATERIAL': (materials.Vacuum, 'Degrader material'),
        'KINETIC_ENERGY': (0.0 * _ureg.MeV, 'Incoming beam energy'),
        'L': (0.0 * _ureg.m, 'Degrader length'),
        'APERTYPE': (None, 'Aperture type (CIRCULAR, ELIPTIC or RECTANGULAR)'),
        'APERTURE': ([], ''),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> List[float]:
        fe = self.material.scattering(kinetic_energy=self.KINETIC_ENERGY,
                                      thickness=self.L,
                                      compute_a1=False,
                                      compute_a2=False)
        return [
            fe['A'][0],
        ]

    def propagate(self,
                  beam_in: _np.ndarray,
                  beam_out: _np.ndarray = None,
                  global_parameters: list = None,
                  ) -> Tuple[_np.ndarray, _np.ndarray]:
        if self.material is materials.Vacuum:
            _np.copyto(dst=beam_out, src=beam_in, casting='no')
            return beam_in, beam_out
        a0 = self.cache[0]

        beam_out[:, 0] = beam_in[:, 0]
        beam_out[:, 2] = beam_in[:, 2]
        beam_out[:, 4] = beam_in[:, 4]
        beam_out[:, 5] = beam_in[:, 5]

        beam_out[:, 1] = beam_in[:, 1] + nprandom.normal(0.0, a0, size=beam_in.shape[0])
        beam_out[:, 3] = beam_in[:, 3] + nprandom.normal(0.0, a0, size=beam_in.shape[0])

        beam_out = self.pseudo_aperture_check(beam_out, self.beta)

        return beam_in, beam_out


class Degrader(MaterialElement):
    PARAMETERS = {
        'MATERIAL': (materials.Vacuum, 'Degrader material'),
        'KINETIC_ENERGY': (0.0 * _ureg.MeV, 'Incoming beam energy'),
        'L': (0.0 * _ureg.m, 'Degrader length'),
        'WITH_LOSSES': (False, 'Boolean to compute losses and dpp'),
        'APERTYPE': (None, 'Aperture type (CIRCULAR, ELIPTIC or RECTANGULAR)'),
        'APERTURE': ([], ''),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> List[float]:
        fe = self.material.scattering(kinetic_energy=self.KINETIC_ENERGY, thickness=self.L)
        return [
            self.L.m_as('m'),
            fe['A'][0],
            fe['A'][1],
            fe['A'][2],
            self.material.energy_dispersion(energy=self.degraded_energy),
            self.material.losses(energy=self.degraded_energy)
        ]

    def propagate(self,
                  beam_in: _np.ndarray,
                  beam_out: _np.ndarray = None,
                  global_parameters: list = None,
                  ) -> Tuple[_np.ndarray, _np.ndarray]:
        length, a0, a1, a2, dpp, losses = self.cache

        if length == 0:
            _np.copyto(dst=beam_out, src=beam_in, casting='no')
            return beam_in, beam_out

        # Monte-Carlo method
        # Remove particles
        if self.WITH_LOSSES is True and losses != 1:
            idx = _np.random.randint(beam_in.shape[0], size=int(losses * beam_in.shape[0]))
            beam_out = beam_in[idx, :]
        else:
            beam_out = beam_in

        # Transport matrix
        matrix = _np.array(
            [
                [1, length, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, length, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
            ]
        )
        beam_out = beam_out.dot(matrix.T)

        # If the integrator is not MAD-X based, swap the last 2 columns
        if self.integrator not in [MadXIntegrator, MadXParaxialDriftIntegrator]:
            beam_out[:, [-2, -1]] = beam_out[:, [-1, -2]]

        # Interactions
        if self.material is not materials.Vacuum:
            beam_out += nprandom.multivariate_normal(
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                _np.array(
                    [
                        [a2, a1, 0.0, 0.0, 0.0, 0.0],
                        [a1, a0, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, a2, a1, 0.0, 0.0],
                        [0.0, 0.0, a1, a0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, dpp ** 2, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ]
                ),
                int(beam_out.shape[0]))

        if self.integrator not in [MadXIntegrator, MadXParaxialDriftIntegrator]:
            beam_out[:, [-1, -2]] = beam_out[:, [-2, -1]]
        else:
            beam_out[:, -1] = _Beam.compute_pt(beam_out[:, -2], self.beta, first_order=False)
            beam_out = self.pseudo_aperture_check(beam_out, self.beta)

        return beam_in, beam_out


class BeamStop(MaterialElement):
    PARAMETERS = {
        'MATERIAL': (materials.Vacuum, 'Beam Stop material'),
        'L': (0.0 * _ureg.m, 'Beam Stop length'),
        'RADIUS': (0.0 * _ureg.m, 'Beam Stop radius'),
        'KINETIC_ENERGY': (0.0 * _ureg.MeV, 'Incoming beam energy'),
    }

    @property
    def parameters(self) -> list:
        return [
            self.L.m_as('m'),
            self.RADIUS.m_as('m')
        ]

    def propagate(self,
                  beam_in: _np.ndarray,
                  beam_out: _np.ndarray = None,
                  global_parameters: list = None,
                  ) -> Tuple[_np.ndarray, _np.ndarray]:
        length, radius = self.parameters

        if length == 0 or radius == 0:
            _np.copyto(dst=beam_out, src=beam_in, casting='no')
            return beam_in, beam_out

        else:
            beam_out = _np.compress(
                (beam_in[:, 0] ** 2 + beam_in[:, 2] ** 2) > radius ** 2,
                beam_in,
                axis=0,
            )
            return beam_in, beam_out


class Gap(Degrader):
    pass