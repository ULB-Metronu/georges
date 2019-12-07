"""
TODO
"""
from typing import Tuple
import numpy as _np
try:
    import numpy.random_intel as nprandom
except ModuleNotFoundError:
    import numpy.random as nprandom
from ... import ureg as _ureg
from .elements import ManzoniElement as _ManzoniElement
from ...fermi import materials


class Scatterer(_ManzoniElement):
    pass


class Degrader(_ManzoniElement):
    PARAMETERS = {
        'MATERIAL': (materials.Vacuum, 'Degrader material'),
        'KINETIC_ENERGY': (100 * _ureg.MeV, 'Incoming beam energy'),
        'L': (0.0 * _ureg.m, 'Degrader length'),
    }
    """Parameters of the element, with their default value and their descriptions."""

    @property
    def parameters(self) -> list:
        fe = self.MATERIAL.scattering(kinetic_energy=self.KINETIC_ENERGY, thickness=self.L)
        return [
            self.L.m_as('m'),
            fe['A'][0],
            fe['A'][1],
            fe['A'][2],
        ]

    def propagate(self,
                  beam_in: _np.ndarray,
                  beam_out: _np.ndarray = None,
                  global_parameters: list = None,
                  ) -> Tuple[_np.ndarray, _np.ndarray]:
        length, a0, a1, a2 = self.parameters

        if length == 0:
            _np.copyto(dst=beam_out, src=beam_in, casting='no')
            return beam_in, beam_out

        # Monte-Carlo method
        # Remove particles
        # if e[INDEX_FE_LOSS] != 0:
        #     idx = np.random.randint(b.shape[0], size=int((e[INDEX_FE_LOSS]) * b.shape[0]))
        #     b = b[idx, :]

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
        beam_out = beam_in.dot(matrix.T)

        # Interactions
        beam_out += nprandom.multivariate_normal(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            _np.array(
                [
                    [a2, a1, 0.0, 0.0, 0.0, 0.0],
                    [a1, a0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, a2, a1, 0.0, 0.0],
                    [0.0, 0.0, a1, a0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                ]
            ),
            int(beam_in.shape[0]))
        return beam_in, beam_out
