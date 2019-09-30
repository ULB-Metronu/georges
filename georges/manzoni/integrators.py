from .maps import compute_mad_combined_dipole_matrix, \
    compute_mad_combined_dipole_tensor, \
    compute_mad_quadrupole_matrix, \
    compute_mad_quadrupole_tensor, \
    compute_transport_combined_dipole_matrix, \
    compute_transport_combined_dipole_tensor, \
    compute_transport_quadrupole_matrix, \
    compute_transport_quadrupole_tensor
from .kernels import batched_vector_matrix, batched_vector_matrix_tensor


class IntegratorType(type):
    pass


class Integrator(metaclass=IntegratorType):
    @classmethod
    def propagate(cls, element, beam_in, beam_out):
        pass


class MadIntegrator(Integrator):
    pass


class FirstOrderTaylorMadIntegrator(MadIntegrator):
    MATRICES = {
        'BEND': compute_mad_combined_dipole_matrix,
        'SBEND': compute_mad_combined_dipole_matrix,
        'QUADRUPOLE': compute_mad_quadrupole_matrix
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out):
        return batched_vector_matrix(beam_in, beam_out, *element.cache)

    @classmethod
    def cache(cls, element):
        return [
            cls.MATRICES.get(element.__class__.__name__.upper())(*element.parameters)
        ]


class SecondOrderTaylorMadIntegrator(FirstOrderTaylorMadIntegrator):
    TENSORS = {
        'BEND': compute_mad_combined_dipole_tensor,
        'SBEND': compute_mad_combined_dipole_tensor,
        'QUADRUPOLE': compute_mad_quadrupole_tensor,
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


class TransportIntegrator(Integrator):
    pass


class FirstOrderTransportIntegrator(TransportIntegrator):
    MATRICES = {
        'BEND': compute_transport_combined_dipole_matrix,
        'QUADRUPOLE': compute_transport_quadrupole_matrix,
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
