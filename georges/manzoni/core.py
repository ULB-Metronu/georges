from __future__ import annotations
from typing import TYPE_CHECKING, Optional, List
import numpy as _np
import pandas as _pd
from numba.typed import List as nList
from georges_core.sequences import BetaBlock as _BetaBlock
from georges_core.twiss import \
    compute_alpha_from_matrix, \
    compute_beta_from_matrix, \
    compute_dispersion_from_matrix, \
    compute_dispersion_prime_from_matrix, \
    compute_gamma_from_matrix, \
    compute_jacobian_from_matrix, \
    compute_mu_from_matrix, \
    compute_periodic_twiss

from .beam import Beam as _Beam
from .observers import BeamObserver as _BeamObserver
if TYPE_CHECKING:
    from .input import Input as _Input
    from .observers import Observer as _Observer
    from .. import Kinematics as _Kinematics


def track(beamline: _Input,
          beam: _Beam,
          observers: List[Optional[_Observer]] = None,
          check_apertures_exit: bool = True,
          check_apertures_entry: bool = False
          ):
    """
    Args:
        beamline:
        beam:
        observers:
        check_apertures_exit:
        check_apertures_entry:
    Returns:
    """
    if observers is None:
        observers = []
    global_parameters = nList()
    global_parameters.append(beam.kinematics.beta)
    b1 = _np.copy(beam.distribution)
    b2 = _np.zeros(b1.shape)
    for e in beamline.sequence:
        if check_apertures_entry:
            b2, b1 = e.check_aperture(b2, b1)
            if b1.shape != b2.shape:
                b1 = _np.zeros(b2.shape)
            if b2.shape[0] == 0:
                break
        b1, b2 = e.propagate(b1, b2, global_parameters)
        if check_apertures_exit:
            b1, b2 = e.check_aperture(b1, b2)
        for o in observers:
            if o is not None:
                o(e, b1, b2)
        if b1.shape != b2.shape:
            b1 = _np.zeros(b2.shape)
        if b2.shape[0] == 0:
            break
        b2, b1 = b1, b2


def twiss(beamline: _Input,
          kinematics: _Kinematics,
          reference_particle: _np.ndarray = None,
          offsets=None,
          with_twiss_parametrization: bool = True,
          twiss_init: _BetaBlock = None,
          ) -> _pd.DataFrame:
    """

    Args:
        beamline:
        kinematics:
        reference_particle:
        offsets:
        with_twiss_parametrization:
        twiss_init:

    Returns:

    """
    def track_for_twiss() -> _pd.DataFrame:
        nonlocal reference_particle
        nonlocal offsets

        if reference_particle is None:
            reference_particle = _np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        if offsets is None:
            offsets = _np.array([0.01, 0.01, 0.01, 0.01, 0.01])
        coordinates = _np.array([
            reference_particle,
            reference_particle + [offsets[0], 0.0, 0.0, 0.0, 0.0, 0.0],
            reference_particle + [0.0, offsets[0], 0.0, 0.0, 0.0, 0.0],
            reference_particle + [0.0, 0.0, offsets[0], 0.0, 0.0, 0.0],
            reference_particle + [0.0, 0.0, 0.0, offsets[0], 0.0, 0.0],
            reference_particle + [0.0, 0.0, 0.0, 0.0, 0.0, offsets[0]],
            reference_particle + [-offsets[0], 0.0, 0.0, 0.0, 0.0, 0.0],
            reference_particle + [0.0, -offsets[0], 0.0, 0.0, 0.0, 0.0],
            reference_particle + [0.0, 0.0, -offsets[0], 0.0, 0.0, 0.0],
            reference_particle + [0.0, 0.0, 0.0, -offsets[0], 0.0, 0.0],
            reference_particle + [0.0, 0.0, 0.0, 0.0, 0.0, -offsets[0]],
        ])
        beam = _Beam(kinematics=kinematics, distribution=coordinates)
        observer = _BeamObserver(with_input_beams=False)
        track(beam=beam, beamline=beamline, observers=observer)
        return observer.to_df()

    def compute_matrix_for_twiss(data: _pd.DataFrame) -> _pd.DataFrame:
        normalization = 2 * offsets
        _matrix = {}
        for _, d in data.iterrows():
            label = d['LABEL1']
            m = d['BEAM_OUT']
            m[:, 4] = m[:, 5]
            m = m[:, 0:5]
            _matrix[label] = {
                f'R{i + 1}{j + 1}': (m[i + 1, j] - m[i + 1 + 5, j]) / normalization[i]
                for j in range(0, 5)
                for i in range(0, 5)
            }
        return _pd.DataFrame.from_dict(_matrix, orient='index')

    def compute_parametrization_for_twiss(m: _pd.DataFrame,
                                          with_phase_unrolling: bool = True) -> _pd.DataFrame:
        nonlocal twiss_init
        if twiss_init is None:
            twiss_init = compute_periodic_twiss(m)

        m['BETA11'] = compute_beta_from_matrix(m, twiss_init)
        m['BETA22'] = compute_beta_from_matrix(m, twiss_init, plane=2)
        m['ALPHA11'] = compute_alpha_from_matrix(m, twiss_init)
        m['ALPHA22'] = compute_alpha_from_matrix(m, twiss_init, plane=2)
        m['GAMMA11'] = compute_gamma_from_matrix(m, twiss_init)
        m['GAMMA22'] = compute_gamma_from_matrix(m, twiss_init, plane=2)
        m['MU1'] = compute_mu_from_matrix(m, twiss_init)
        m['MU2'] = compute_mu_from_matrix(m, twiss_init, plane=2)
        m['DET1'] = compute_jacobian_from_matrix(m)
        m['DET2'] = compute_jacobian_from_matrix(m, plane=2)
        m['DISP1'] = compute_dispersion_from_matrix(m, twiss_init)
        m['DISP2'] = compute_dispersion_prime_from_matrix(m, twiss_init)
        m['DISP3'] = compute_dispersion_from_matrix(m, twiss_init, plane=2)
        m['DISP4'] = compute_dispersion_prime_from_matrix(m, twiss_init, plane=2)

        def phase_unrolling(phi):
            """TODO"""
            if phi[0] < 0:
                phi[0] += 2 * _np.pi
            for i in range(1, phi.shape[0]):
                if phi[i] < 0:
                    phi[i] += 2 * _np.pi
                if phi[i - 1] - phi[i] > 0.5:
                    phi[i:] += 2 * _np.pi
            return phi

        try:
            from numba import njit
            phase_unrolling = njit(phase_unrolling)
        except ModuleNotFoundError:
            pass

        if with_phase_unrolling:
            m['MU1U'] = phase_unrolling(m['MU1'].values)
            m['MU2U'] = phase_unrolling(m['MU2'].values)

        return m

    tracks = track_for_twiss()
    matrix = compute_matrix_for_twiss(tracks)
    if with_twiss_parametrization:
        matrix = compute_parametrization_for_twiss(matrix)
    return matrix


def match(beamline: _Input, beam: _Beam):
    ...
