from typing import List
import numpy as _np
from .maps import compute_mad_combined_dipole_matrix, \
    compute_mad_combined_dipole_tensor, \
    compute_mad_quadrupole_matrix, \
    compute_mad_quadrupole_tensor, \
    track_madx_quadrupole, \
    track_madx_drift
from .kernels import batched_vector_matrix, batched_vector_matrix_tensor

__ALL__ = [
    'IntegratorType',
    'Integrator',
    'MadXIntegrator',
    'Mad8Integrator',
    'Mad8FirstOrderTaylorIntegrator',
    'Mad8SecondOrderTaylorIntegrator',
    'TransportIntegrator',
    'TransportFirstOrderIntegrator',
    'TransportSecondOrderIntegrator',
    'PTCIntegrator',
]


class IntegratorType(type):
    pass


class Integrator(metaclass=IntegratorType):
    @classmethod
    def propagate(cls, element, beam_in: _np.ndarray, beam_out: _np.ndarray, global_parameters: list):
        return beam_in, beam_out


class MadXIntegrator(Integrator):
    METHODS = {
        'BEND': None,
        'SBEND': None,
        'DRIFT': track_madx_drift,
        'QUADRUPOLE': track_madx_quadrupole,
    }

    @classmethod
    def propagate(cls, element, beam_in: _np.ndarray, beam_out: _np.ndarray, global_parameters: list):
        return cls.METHODS.get(element.__class__.__name__.upper())(
            beam_in, beam_out, element.cache, global_parameters
        )

    @classmethod
    def cache(cls, element) -> list:
        return element.parameters


class Mad8Integrator(Integrator):
    pass


class Mad8FirstOrderTaylorIntegrator(Mad8Integrator):
    MATRICES = {
        'BEND': compute_mad_combined_dipole_matrix,
        'SBEND': compute_mad_combined_dipole_matrix,
        'QUADRUPOLE': compute_mad_quadrupole_matrix,
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out):
        return batched_vector_matrix(beam_in, beam_out, *element.cache)

    @classmethod
    def cache(cls, element) -> List:
        return [
            cls.MATRICES.get(element.__class__.__name__.upper())(*element.parameters)
        ]


class Mad8SecondOrderTaylorIntegrator(Mad8FirstOrderTaylorIntegrator):
    TENSORS = {
        'BEND': compute_mad_combined_dipole_tensor,
        'SBEND': compute_mad_combined_dipole_tensor,
        'QUADRUPOLE': compute_mad_quadrupole_tensor,
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out):
        return batched_vector_matrix_tensor(beam_in, beam_out, *element.cache)

    @classmethod
    def cache(cls, element) -> List:
        return [
            cls.MATRICES.get(element.__class__.__name__.upper())(*element.parameters),
            cls.TENSORS.get(element.__class__.__name__.upper())(*element.parameters)
        ]


class TransportIntegrator(Integrator):
    pass


class TransportFirstOrderIntegrator(TransportIntegrator):
    pass


class TransportSecondOrderIntegrator(TransportFirstOrderIntegrator):
    pass


class PTCIntegrator(Integrator):
    pass
