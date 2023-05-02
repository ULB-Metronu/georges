"""
The file `integrators.py` contains the definition of all the integrator types,
with a propagate function that is a wrapper to the `propagate` function of each of the
different elements in the beamline. Each integrator has a dictionary with the different elements
it can be used with, so that if this integrator is selected, the propagate function of the element
will be called via the propagate function of the integrator to allow the tracking.

"""
from typing import List, Tuple

import numpy as _np
from numba.typed import List as nList

from .kernels import batched_vector_matrix, batched_vector_matrix_tensor
from .maps import (
    compute_mad_combined_dipole_matrix,
    compute_mad_combined_dipole_tensor,
    compute_mad_drift_matrix,
    compute_mad_drift_tensor,
    compute_mad_quadrupole_matrix,
    compute_mad_quadrupole_tensor,
    compute_transport_combined_dipole_ex_matrix,
    compute_transport_combined_dipole_ex_tensor,
    compute_transport_combined_dipole_matrix,
    compute_transport_combined_dipole_tensor,
    compute_transport_fringe_in_ex_matrix,
    compute_transport_fringe_in_ex_tensor,
    compute_transport_fringe_in_matrix,
    compute_transport_fringe_in_tensor,
    compute_transport_fringe_out_ex_matrix,
    compute_transport_fringe_out_ex_tensor,
    compute_transport_fringe_out_matrix,
    compute_transport_fringe_out_tensor,
    compute_transport_multipole_ex_matrix,
    compute_transport_multipole_ex_tensor,
    compute_transport_multipole_matrix,
    compute_transport_multipole_tensor,
    compute_transport_quadrupole_ex_matrix,
    compute_transport_quadrupole_ex_tensor,
    compute_transport_quadrupole_matrix,
    compute_transport_quadrupole_tensor,
    compute_transport_sextupole_ex_matrix,
    compute_transport_sextupole_ex_tensor,
    compute_transport_sextupole_matrix,
    compute_transport_sextupole_tensor,
    track_madx_bend,
    track_madx_dipedge,
    track_madx_drift,
    track_madx_drift_paraxial,
    track_madx_kicker,
    track_madx_quadrupole,
    track_madx_srotation,
)

__all__ = [
    "IntegratorType",
    "Integrator",
    "MadXIntegrator",
    "MadXParaxialDriftIntegrator",
    "Mad8Integrator",
    "Mad8FirstOrderTaylorIntegrator",
    "Mad8SecondOrderTaylorIntegrator",
    "TransportIntegrator",
    "TransportFirstOrderTaylorIntegrator",
    "TransportFirstOrderTaylorIntegratorExact",
    "TransportSecondOrderTaylorIntegrator",
    "TransportSecondOrderTaylorIntegratorExact",
    "PTCIntegrator",
]


class IntegratorType(type):
    pass


class Integrator(metaclass=IntegratorType):
    @classmethod
    def propagate(
        cls,
        element,
        beam_in: _np.ndarray,
        beam_out: _np.ndarray,
        global_parameters: nList,
    ) -> Tuple[_np.ndarray, _np.ndarray]:
        return beam_in, beam_out


class MadXIntegrator(Integrator):
    METHODS = {
        "DIPEDGE": track_madx_dipedge,
        "RBEND": track_madx_bend,
        "SBEND": track_madx_bend,
        "DRIFT": track_madx_drift,
        "GAP": track_madx_drift,
        "RECTANGULARCOLLIMATOR": track_madx_drift,
        "ELLIPTICALCOLLIMATOR": track_madx_drift,
        "CIRCULARCOLLIMATOR": track_madx_drift,
        "DUMP": track_madx_drift,
        "QUADRUPOLE": track_madx_quadrupole,
        "SROTATION": track_madx_srotation,
        "KICKER": track_madx_kicker,
        "HKICKER": track_madx_kicker,
        "VKICKER": track_madx_kicker,
    }

    @classmethod
    def propagate(cls, element, beam_in: _np.ndarray, beam_out: _np.ndarray, global_parameters: nList):
        return cls.METHODS.get(element.__class__.__name__.upper())(
            beam_in,
            beam_out,
            element.cache,
            global_parameters,
        )

    @classmethod
    def cache(cls, element) -> List:
        return element.parameters


class MadXParaxialDriftIntegrator(MadXIntegrator):
    METHODS = {
        "DIPEDGE": track_madx_dipedge,
        "RBEND": track_madx_bend,
        "SBEND": track_madx_bend,
        "DRIFT": track_madx_drift_paraxial,
        "RECTANGULARCOLLIMATOR": track_madx_drift_paraxial,
        "ELLIPTICALCOLLIMATOR": track_madx_drift_paraxial,
        "CIRCULARCOLLIMATOR": track_madx_drift_paraxial,
        "DUMP": track_madx_drift_paraxial,
        "QUADRUPOLE": track_madx_quadrupole,
        "SROTATION": track_madx_srotation,
        "KICKER": track_madx_kicker,
        "HKICKER": track_madx_kicker,
        "VKICKER": track_madx_kicker,
    }


class Mad8Integrator(Integrator):
    pass


class Mad8FirstOrderTaylorIntegrator(Mad8Integrator):
    MATRICES = {
        "DRIFT": compute_mad_drift_matrix,
        "BEND": compute_mad_combined_dipole_matrix,
        "SBEND": compute_mad_combined_dipole_matrix,
        "QUADRUPOLE": compute_mad_quadrupole_matrix,
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out, global_parameters: nList):
        if element.__class__.__name__.upper() in ["HKICKER", "VKICKER"]:
            return track_madx_kicker(beam_in, beam_out, element.cache, global_parameters)

        elif element.__class__.__name__.upper() in [
            "GAP",
            "RECTANGULARCOLLIMATOR",
            "ELLIPTICALCOLLIMATOR",
            "CIRCULARCOLLIMATOR",
            "DUMP",
        ]:
            return batched_vector_matrix(
                beam_in,
                beam_out,
                compute_mad_drift_matrix(element.cache, global_parameters),
            )

        elif element.__class__.__name__.upper() == "SROTATION":
            return track_madx_srotation(beam_in, beam_out, element.cache, global_parameters)

        else:
            return batched_vector_matrix(
                beam_in,
                beam_out,
                cls.MATRICES.get(element.__class__.__name__.upper())(element.cache, global_parameters),
            )

    @classmethod
    def cache(cls, element) -> List:
        return element.parameters


class Mad8SecondOrderTaylorIntegrator(Mad8FirstOrderTaylorIntegrator):
    TENSORS = {
        "DRIFT": compute_mad_drift_tensor,
        "BEND": compute_mad_combined_dipole_tensor,
        "SBEND": compute_mad_combined_dipole_tensor,
        "QUADRUPOLE": compute_mad_quadrupole_tensor,
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out, global_parameters: nList):
        if element.__class__.__name__.upper() in ["HKICKER", "VKICKER"]:
            return track_madx_kicker(beam_in, beam_out, element.cache, global_parameters)

        elif element.__class__.__name__.upper() in [
            "GAP",
            "RECTANGULARCOLLIMATOR",
            "ELLIPTICALCOLLIMATOR",
            "CIRCULARCOLLIMATOR",
            "DUMP",
        ]:
            return batched_vector_matrix_tensor(
                beam_in,
                beam_out,
                compute_mad_drift_matrix(element.cache, global_parameters),
                compute_mad_drift_tensor(element.cache, global_parameters),
            )

        elif element.__class__.__name__.upper() == "SROTATION":
            return track_madx_srotation(beam_in, beam_out, element.cache, global_parameters)

        else:
            return batched_vector_matrix_tensor(
                beam_in,
                beam_out,
                cls.MATRICES.get(element.__class__.__name__.upper())(element.cache, global_parameters),
                cls.TENSORS.get(element.__class__.__name__.upper())(element.cache, global_parameters),
            )

    @classmethod
    def cache(cls, element) -> List:
        return element.parameters


class TransportIntegrator(Integrator):
    pass


class TransportFirstOrderTaylorIntegrator(TransportIntegrator):
    MATRICES = {
        "BEND": compute_transport_combined_dipole_matrix,
        "SBEND": compute_transport_combined_dipole_matrix,
        "QUADRUPOLE": compute_transport_quadrupole_matrix,
        "SEXTUPOLE": compute_transport_sextupole_matrix,
        "MULTIPOLE": compute_transport_multipole_matrix,
        "FRINGEIN": compute_transport_fringe_in_matrix,
        "FRINGEOUT": compute_transport_fringe_out_matrix,
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out, global_parameters: nList):
        if element.__class__.__name__.upper() in ["HKICKER", "VKICKER"]:
            return track_madx_kicker(beam_in, beam_out, element.cache, global_parameters)

        elif element.__class__.__name__.upper() in [
            "DRIFT",
            "GAP",
            "RECTANGULARCOLLIMATOR",
            "ELLIPTICALCOLLIMATOR",
            "CIRCULARCOLLIMATOR",
            "DUMP",
        ]:
            return track_madx_drift(beam_in, beam_out, element.cache, global_parameters)

        elif element.__class__.__name__.upper() == "SROTATION":
            return track_madx_srotation(beam_in, beam_out, element.cache, global_parameters)

        else:
            b = batched_vector_matrix(
                beam_in,
                beam_out,
                cls.MATRICES.get(element.__class__.__name__.upper())(element.cache),
            )
            return b[0], b[1]

    @classmethod
    def cache(cls, element) -> List:
        return element.parameters


class TransportFirstOrderTaylorIntegratorExact(TransportIntegrator):
    MATRICES = {
        "BEND": compute_transport_combined_dipole_ex_matrix,
        "SBEND": compute_transport_combined_dipole_ex_matrix,
        "QUADRUPOLE": compute_transport_quadrupole_ex_matrix,
        "SEXTUPOLE": compute_transport_sextupole_ex_matrix,
        "MULTIPOLE": compute_transport_multipole_ex_matrix,
        "FRINGEIN": compute_transport_fringe_in_ex_matrix,
        "FRINGEOUT": compute_transport_fringe_out_ex_matrix,
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out, global_parameters: nList):
        if element.__class__.__name__.upper() in ["HKICKER", "VKICKER"]:
            return track_madx_kicker(beam_in, beam_out, element.cache, global_parameters)

        elif element.__class__.__name__.upper() in [
            "DRIFT",
            "GAP",
            "RECTANGULARCOLLIMATOR",
            "ELLIPTICALCOLLIMATOR",
            "CIRCULARCOLLIMATOR",
            "DUMP",
        ]:
            return track_madx_drift(beam_in, beam_out, element.cache, global_parameters)

        elif element.__class__.__name__.upper() == "SROTATION":
            return track_madx_srotation(beam_in, beam_out, element.cache, global_parameters)

        else:
            b2 = _np.zeros(beam_in.shape)
            for i in range(beam_in.shape[0]):
                updated_parameters = element.cache.copy()

                updated_parameters.append(beam_in[i, 4])
                matrix = cls.MATRICES.get(element.__class__.__name__.upper())(updated_parameters)
                b2[i, :] = batched_vector_matrix(_np.array([beam_in[i, :]]), _np.array([beam_out[i, :]]), matrix)[1]

            return beam_in, b2

    @classmethod
    def cache(cls, element) -> List:
        return element.parameters


class TransportSecondOrderTaylorIntegrator(TransportFirstOrderTaylorIntegrator):
    TENSORS = {
        "BEND": compute_transport_combined_dipole_tensor,
        "SBEND": compute_transport_combined_dipole_tensor,
        "QUADRUPOLE": compute_transport_quadrupole_tensor,
        "SEXTUPOLE": compute_transport_sextupole_tensor,
        "MULTIPOLE": compute_transport_multipole_tensor,
        "FRINGEIN": compute_transport_fringe_in_tensor,
        "FRINGEOUT": compute_transport_fringe_out_tensor,
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out, global_parameters: nList):
        if element.__class__.__name__.upper() in ["HKICKER", "VKICKER"]:
            return track_madx_kicker(beam_in, beam_out, element.cache, global_parameters)

        elif element.__class__.__name__.upper() in [
            "DRIFT",
            "GAP",
            "RECTANGULARCOLLIMATOR",
            "ELLIPTICALCOLLIMATOR",
            "CIRCULARCOLLIMATOR",
            "DUMP",
        ]:
            return track_madx_drift(beam_in, beam_out, element.cache, global_parameters)

        elif element.__class__.__name__.upper() == "SROTATION":
            return track_madx_srotation(beam_in, beam_out, element.cache, global_parameters)

        else:
            b = batched_vector_matrix_tensor(
                beam_in,
                beam_out,
                cls.MATRICES.get(element.__class__.__name__.upper())(element.cache),
                cls.TENSORS.get(element.__class__.__name__.upper())(element.cache),
            )
            return b[0], b[1]

    @classmethod
    def cache(cls, element) -> List:
        return element.parameters


class TransportSecondOrderTaylorIntegratorExact(TransportFirstOrderTaylorIntegratorExact):
    TENSORS = {
        "BEND": compute_transport_combined_dipole_ex_tensor,
        "SBEND": compute_transport_combined_dipole_ex_tensor,
        "QUADRUPOLE": compute_transport_quadrupole_ex_tensor,
        "SEXTUPOLE": compute_transport_sextupole_ex_tensor,
        "MULTIPOLE": compute_transport_multipole_ex_tensor,
        "FRINGEIN": compute_transport_fringe_in_ex_tensor,
        "FRINGEOUT": compute_transport_fringe_out_ex_tensor,
    }

    @classmethod
    def propagate(cls, element, beam_in, beam_out, global_parameters: nList):
        if element.__class__.__name__.upper() in ["HKICKER", "VKICKER"]:
            return track_madx_kicker(beam_in, beam_out, element.cache, global_parameters)

        elif element.__class__.__name__.upper() in [
            "DRIFT",
            "GAP",
            "RECTANGULARCOLLIMATOR",
            "ELLIPTICALCOLLIMATOR",
            "CIRCULARCOLLIMATOR",
            "DUMP",
        ]:
            return track_madx_drift(beam_in, beam_out, element.cache, global_parameters)

        elif element.__class__.__name__.upper() == "SROTATION":
            return track_madx_srotation(beam_in, beam_out, element.cache, global_parameters)

        else:
            b2 = _np.zeros(beam_in.shape)
            for i in range(beam_in.shape[0]):
                updated_parameters = element.cache.copy()
                updated_parameters.append(beam_in[i, 5])
                matrix = cls.MATRICES.get(element.__class__.__name__.upper())(updated_parameters)
                tensor = cls.TENSORS.get(element.__class__.__name__.upper())(updated_parameters)
                b2[i, :] = batched_vector_matrix_tensor(
                    _np.array([beam_in[i, :]]),
                    _np.array([beam_out[i, :]]),
                    matrix,
                    tensor,
                )[1]

            return beam_in, b2

    @classmethod
    def cache(cls, element) -> List:
        return element.parameters


class PTCIntegrator(Integrator):
    pass
