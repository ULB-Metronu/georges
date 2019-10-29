from typing import List, Tuple
import numpy as _np
from .maps import compute_mad_combined_dipole_matrix, \
    compute_mad_combined_dipole_tensor, \
    compute_mad_quadrupole_matrix, \
    compute_mad_quadrupole_tensor, \
    compute_transport_combined_dipole_matrix, \
    compute_transport_combined_dipole_tensor, \
    compute_transport_multipole_matrix, \
    compute_transport_multipole_tensor, \
    compute_transport_quadrupole_matrix, \
    compute_transport_quadrupole_tensor, \
    compute_transport_sextupole_matrix, \
    compute_transport_sextupole_tensor, \
    track_madx_quadrupole, \
    track_madx_drift, \
    track_madx_bend, \
    track_madx_dipedge
from .kernels import batched_vector_matrix, batched_vector_matrix_tensor, batched_vector_tensor

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
    def propagate(cls,
                  element,
                  beam_in: _np.ndarray,
                  beam_out: _np.ndarray,
                  global_parameters: list
                  ) -> Tuple[_np.ndarray, _np.ndarray]:
        return beam_in, beam_out


class MadXIntegrator(Integrator):
    METHODS = {
        'DIPEDGE': track_madx_dipedge,
        'RBEND': track_madx_bend,
        'SBEND': track_madx_bend,
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
    def propagate(cls, element, beam_in, beam_out, global_parameters: list):
        return batched_vector_matrix(
            beam_in,
            beam_out,
            cls.MATRICES.get(element.__class__.__name__.upper())(*element.cache, *global_parameters)
        )

    @classmethod
    def cache(cls, element) -> List:
        return element.parameters


class Mad8SecondOrderTaylorIntegrator(Mad8FirstOrderTaylorIntegrator):
    TENSORS = {
        'BEND': compute_mad_combined_dipole_tensor,
        'SBEND': compute_mad_combined_dipole_tensor,
        'QUADRUPOLE': compute_mad_quadrupole_tensor,
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out, global_parameters: list):
        return batched_vector_matrix_tensor(
            beam_in,
            beam_out,
            cls.MATRICES.get(element.__class__.__name__.upper())(*element.cache, *global_parameters),
            cls.TENSORS.get(element.__class__.__name__.upper())(*element.cache, *global_parameters)
        )

    @classmethod
    def cache(cls, element) -> List:
        return element.parameters


class TransportIntegrator(Integrator):
    pass


class FirstOrderTransportIntegrator(TransportIntegrator):
    MATRICES = {
        'BEND': compute_transport_combined_dipole_matrix,
        'QUADRUPOLE': compute_transport_quadrupole_matrix,
        'SEXTUPOLE': compute_transport_sextupole_matrix,
        'MULTIPOLE': compute_transport_multipole_matrix,
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out):
        return batched_vector_matrix(beam_in, beam_out, *element.cache)

    @classmethod
    def cache(cls, element):
        return [
            cls.MATRICES.get(element.__class__.__name__.upper())(*element.parameters)
        ]


class SecondOrderTransportIntegrator(FirstOrderTransportIntegrator):
    TENSORS = {
        'BEND': compute_transport_combined_dipole_tensor,
        'QUADRUPOLE': compute_transport_quadrupole_tensor,
        'SEXTUPOLE': compute_transport_sextupole_tensor,
        'MULTIPOLE': compute_transport_multipole_tensor,
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out):
        return batched_vector_matrix_tensor(beam_in, beam_out, *element.cache)

    @classmethod
    def cache(cls, element):
        return [
            cls.MATRICES.get(element.__class__.__name__.upper())(*element.parameters),
            cls.TENSORS.get(element.__class__.__name__.upper())(*element.parameters)
        ]


class PTCIntegrator(Integrator):
    pass
