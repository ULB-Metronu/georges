"""
The `maps` submodule contains the physics of the propagation of the particles through each of the
beamline elements described in the `elements` submodule. For each element, a first and second-order
type propagation is implemented, allowing the user to select the order of the tracking that is suitable
for his specific application. It is done via the selection of the integrator type when building the
beamline (see later). The `Transport-type`, `MadX-type`, and `Mad8-type` maps are implemented and
available. The user should be aware that the canonical variables of the particles are not the same
for the three different integrator types, so the definition of the beam must be done according to
the integrator to be consistent.
"""

from .mad8_combined_dipole import compute_mad_combined_dipole_matrix, compute_mad_combined_dipole_tensor
from .mad8_drift import compute_mad_drift_matrix, compute_mad_drift_tensor, drift6
from .mad8_quadrupole import compute_mad_quadrupole_matrix, compute_mad_quadrupole_tensor
from .madx_combined_dipole import tmsect
from .madx_thick import (
    track_madx_bend,
    track_madx_dipedge,
    track_madx_drift,
    track_madx_drift_paraxial,
    track_madx_kicker,
    track_madx_quadrupole,
    track_madx_srotation,
)
from .transport_combined_dipole import (
    compute_transport_combined_dipole_matrix,
    compute_transport_combined_dipole_tensor,
)
from .transport_combined_dipole_ex import (
    compute_transport_combined_dipole_ex_matrix,
    compute_transport_combined_dipole_ex_tensor,
)
from .transport_drift import compute_transport_drift_matrix
from .transport_fringe_in import compute_transport_fringe_in_matrix, compute_transport_fringe_in_tensor
from .transport_fringe_in_ex import compute_transport_fringe_in_ex_matrix, compute_transport_fringe_in_ex_tensor
from .transport_fringe_out import compute_transport_fringe_out_matrix, compute_transport_fringe_out_tensor
from .transport_fringe_out_ex import compute_transport_fringe_out_ex_matrix, compute_transport_fringe_out_ex_tensor
from .transport_multipole import compute_transport_multipole_matrix, compute_transport_multipole_tensor
from .transport_multipole_ex import compute_transport_multipole_ex_matrix, compute_transport_multipole_ex_tensor
from .transport_quadrupole import compute_transport_quadrupole_matrix, compute_transport_quadrupole_tensor
from .transport_quadrupole_ex import compute_transport_quadrupole_ex_matrix, compute_transport_quadrupole_ex_tensor
from .transport_sextupole import compute_transport_sextupole_matrix, compute_transport_sextupole_tensor
from .transport_sextupole_ex import compute_transport_sextupole_ex_matrix, compute_transport_sextupole_ex_tensor
