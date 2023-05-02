"""
The beam propagation from one element to another in the beamline is implemented in the file `core.py`.
The `track` function contains the loop over the different elements, with optionally the check of
the apertures to select the particles that survive the tracking at the end of each element,
and the use of “observers” if defined by the user, to save the data during the tracking.
The file `core.py` also contains a Twiss function, which allows the calculation of the matrix elements
of all the elements along a given beamline for the Twiss functions calculation based on the 11 particles
method. The user must be aware that this function also needs the georges_core module to work properly,
as the Twiss computation is done in the end in the georges_core library, using the matrix elements
calculated in georges.
    """

from __future__ import annotations

from typing import TYPE_CHECKING, List, Optional

import numpy as _np
import pandas as _pd
from georges_core.sequences import BetaBlock as _BetaBlock
from georges_core.twiss import Twiss as _Twiss
from numba.typed import List as nList

from .beam import Beam as _Beam
from .observers import BeamObserver as _BeamObserver

if TYPE_CHECKING:
    from .. import Kinematics as _Kinematics
    from .input import Input as _Input
    from .observers import Observer as _Observer


def track(
    beamline: _Input,
    beam: _Beam,
    observers: List[Optional[_Observer]] = None,
    check_apertures_exit: bool = False,
    check_apertures_entry: bool = False,
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


def twiss(
    beamline: _Input,
    kinematics: _Kinematics,
    reference_particle: _np.ndarray = None,
    offsets=None,
    twiss_parametrization: bool = True,
    twiss_init: _BetaBlock = None,
    with_phase_unrolling: bool = True,
) -> _pd.DataFrame:
    """

    Args:
        beamline:
        kinematics:
        reference_particle:
        offsets:
        twiss_parametrization:
        twiss_init:
        with_phase_unrolling:
    Returns:

    """

    def track_for_twiss() -> _pd.DataFrame:
        nonlocal reference_particle
        nonlocal offsets

        if reference_particle is None:
            reference_particle = _np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        if offsets is None:
            offsets = _np.array([0.01, 0.01, 0.01, 0.01, 0.01])
        pt = _Beam.compute_pt(dpp=offsets[4], beta=kinematics.beta)
        coordinates = _np.array(
            [
                reference_particle,
                reference_particle + [offsets[0], 0.0, 0.0, 0.0, 0.0, 0.0],
                reference_particle + [0.0, offsets[1], 0.0, 0.0, 0.0, 0.0],
                reference_particle + [0.0, 0.0, offsets[2], 0.0, 0.0, 0.0],
                reference_particle + [0.0, 0.0, 0.0, offsets[3], 0.0, 0.0],
                reference_particle + [0.0, 0.0, 0.0, 0.0, offsets[4], pt],
                reference_particle + [-offsets[0], 0.0, 0.0, 0.0, 0.0, 0.0],
                reference_particle + [0.0, -offsets[1], 0.0, 0.0, 0.0, 0.0],
                reference_particle + [0.0, 0.0, -offsets[2], 0.0, 0.0, 0.0],
                reference_particle + [0.0, 0.0, 0.0, -offsets[3], 0.0, 0.0],
                reference_particle + [0.0, 0.0, 0.0, 0.0, -offsets[4], -pt],
            ],
        )
        beam = _Beam(kinematics=kinematics, distribution=coordinates)
        observer = _BeamObserver(with_input_beams=False)
        track(beam=beam, beamline=beamline, observers=[observer])
        return observer.to_df()

    def compute_matrix_for_twiss(data: _pd.DataFrame) -> _pd.DataFrame:
        # Make some changes in the format of the matrices to be consistent with georges-core
        # Indeed in georges core, we have a longitudinal coordinates
        normalization = 2 * offsets
        normalization = _np.insert(normalization, -1, normalization[0], 0)
        _matrix = {}
        for label, d in data.iterrows():
            m = d["BEAM_OUT"]
            m[:, 4] = m[:, 5]
            m = m[:, 0:5]
            m = _np.hstack((m, _np.zeros((m.shape[0], 1))))
            m[:, [-2, -1]] = m[:, [-1, -2]]
            m = _np.insert(m, 5, _np.zeros(m.shape[1]), 0)
            m = _np.insert(m, -1, _np.zeros(m.shape[1]), 0)
            _matrix[label] = {
                f"R{j + 1}{i + 1}": (m[i + 1, j] - m[i + 1 + 6, j]) / normalization[i]
                for j in range(0, 6)
                for i in range(0, 6)
            }
            _matrix[label]["S"] = d["AT_EXIT"].m_as("m")
        return _pd.DataFrame.from_dict(_matrix, orient="index")

    tracks = track_for_twiss()
    matrix = compute_matrix_for_twiss(tracks)
    if twiss_parametrization:
        matrix = _Twiss(twiss_init=twiss_init, with_phase_unrolling=with_phase_unrolling)(matrix=matrix)
    return matrix


def match(beamline: _Input, beam: _Beam):
    ...
