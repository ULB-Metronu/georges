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


class Scatterer(_ManzoniElement):
    PARAMETERS = {
        'MATERIAL': (materials.Vacuum, 'Degrader material'),
        'KINETIC_ENERGY': (82.5 * _ureg.MeV, 'Incoming beam energy'),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> List[float]:
        fe = self.MATERIAL.scattering(kinetic_energy=self.KINETIC_ENERGY, thickness=self.L)
        return [
            fe['A'][0],
        ]

    def propagate(self,
                  beam_in: _np.ndarray,
                  beam_out: _np.ndarray = None,
                  global_parameters: list = None,
                  ) -> Tuple[_np.ndarray, _np.ndarray]:
        a0 = self.parameters[0]

        beam_out[:, 1] = beam_in[:, 1] + nprandom.normal(0.0, a0, size=beam_in.shape[0])
        beam_out[:, 3] = beam_in[:, 3] + nprandom.normal(0.0, a0, size=beam_in.shape[0])

        return beam_in, beam_out


class Degrader(_ManzoniElement):
    PARAMETERS = {
        'MATERIAL': (materials.Vacuum, 'Degrader material'),
        'KINETIC_ENERGY': (82.5 * _ureg.MeV, 'Incoming beam energy'),
        'L': (0.0 * _ureg.m, 'Degrader length'),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> List[float]:
        fe = self.MATERIAL.scattering(kinetic_energy=self.KINETIC_ENERGY, thickness=self.L)
        return [
            self.L.m_as('m'),
            fe['A'][0],
            fe['A'][1],
            fe['A'][2],
            self.MATERIAL.energy_dispersion(energy=self.KINETIC_ENERGY),
            self.MATERIAL.losses(energy=self.KINETIC_ENERGY)
        ]

    def propagate(self,
                  beam_in: _np.ndarray,
                  beam_out: _np.ndarray = None,
                  global_parameters: list = None,
                  ) -> Tuple[_np.ndarray, _np.ndarray]:
        length, a0, a1, a2, dpp, losses = self.parameters

        if length == 0:
            _np.copyto(dst=beam_out, src=beam_in, casting='no')
            return beam_in, beam_out

        # Monte-Carlo method
        # Remove particles
        if losses != 0:
            idx = _np.random.randint(beam_in.shape[0], size=int(losses * beam_in.shape[0]))
            beam_out = beam_in[idx, :]

        # Transport
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

        # Interactions
        beam_out += nprandom.multivariate_normal(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            _np.array(
                [
                    [a2, a1, 0.0, 0.0, 0.0, 0.0],
                    [a1, a0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, a2, a1, 0.0, 0.0],
                    [0.0, 0.0, a1, a0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, dpp**2, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                ]
            ),
            int(beam_out.shape[0]))

        return beam_in, beam_out


class BeamStop(_ManzoniElement):
    PARAMETERS = {
        'MATERIAL': (materials.Lead, 'Beam Stop material'),
        'L': (0.0 * _ureg.m, 'Beam Stop length'),
        'RADIUS': (0.0 * _ureg.m, 'Beam Stop radius'),
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
