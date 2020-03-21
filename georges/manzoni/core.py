from __future__ import annotations
from typing import TYPE_CHECKING, Optional
import numpy as _np
import pandas as _pd
if TYPE_CHECKING:
    from .input import Input as _Input
    from .beam import Beam as _Beam
    from .observers import Observer as _Observer
    from .observers import BeamObserver as _BeamObserver
    from .. import Kinematics as _Kinematics


def track(beamline: _Input,
          beam: _Beam,
          observer: Optional[_Observer] = None,
          check_apertures_exit: bool = True,
          check_apertures_entry: bool = False
          ):
    """
    Args:
        beamline:
        beam:
        observer:
        check_apertures_exit:
        check_apertures_entry:
    Returns:
    """
    global_parameters = [
        beam.kinematics.beta
    ]
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
        if observer is not None:
            observer(e, b1, b2)
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
          twiss_init: _pd.Series = None,
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
        track(beam=beam, beamline=beamline, observer=observer)
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

    def compute_parametrization_for_twiss(m: _pd.DataFrame) -> _pd.DataFrame:
        def compute_periodic_twiss() -> _pd.Series:
            nonlocal m
            # Here compute periodic twiss based on matrix m
            return _pd.Series()

        nonlocal twiss_init
        if twiss_init is None:
            twiss_init = compute_periodic_twiss()

        # Here compute twiss parametrization based on twiss_init and matrix m

        return m

    tracks = track_for_twiss()
    matrix = compute_matrix_for_twiss(tracks)
    if with_twiss_parametrization:
        matrix = compute_parametrization_for_twiss(matrix)
    return matrix


def match(beamline: _Input, beam: _Beam):
    ...
