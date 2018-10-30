import numpy as np
from .matrices import matrices, matrices4
from .tensors import tensors
from .integrators import integrators
from .fe import mc, fe
from .constants import *
from .aperture import aperture_check


class ManzoniException(Exception):
    """Exception raised for errors in the Manzoni module."""

    def __init__(self, m):
        self.message = m


def sigma(line, beam, **kwargs):
    """
    Sigma-matrix tracking through a beamline.
    Code optimized for performance.
    :param line: beamline description in Manzoni format
    :param beam: initial beam sigma matrix
    :param kwargs: optional parameters
    :return: Observer.track_end() return value
    """
    sigmas = np.zeros([line.shape[0], 4, 4])
    sigmas[0, :] = beam
    # Main loop
    for i in range(0, line.shape[0]):
        # Track
        if line[i, INDEX_CLASS_CODE] in CLASS_CODE_MATRIX:
            matrix = matrices4[int(line[i, INDEX_CLASS_CODE])](line[i])
            beam = np.matmul(np.matmul(matrix, beam), matrix.T)
        elif line[i, INDEX_CLASS_CODE] in CLASS_CODE_FE:
            beam = fe[int(line[i, INDEX_CLASS_CODE])](line[i], beam, **kwargs)
        else:
            # Default to a drift
            matrix = matrices4[CLASS_CODES['DRIFT']](line[i])
            beam = np.matmul(np.matmul(matrix, beam), matrix.T)
        sigmas[i, :] = beam
    return sigmas


def track(line, beam, turns=1, observer=None, order=1, **kwargs):
    """
    Tracking through a beamline.
    Code optimized for performance.
    :param line: beamline description in Manzoni format
    :param beam: initial beam
    :param turns: number of tracking turns
    :param observer: Observer object to witness and record the tracking data
    :param order: Integration order (default: 1)
    :param kwargs: optional parameters
    :return: Observer.track_end() return value
    """
    if order == 1:
        return track1(line, beam, turns, observer, **kwargs)
    elif order == 2:
        return track2(line, beam, turns, observer, **kwargs)
    else:
        raise ManzoniException("Invalid tracking order")


def track1(line, beam, turns=1, observer=None, **kwargs):
    """
    Tracking through a beamline.
    Code optimized for performance.
    :param line: beamline description in Manzoni format
    :param beam: initial beam
    :param turns: number of tracking turns
    :param observer: Observer object to witness and record the tracking data
    :param kwargs: optional parameters
    :return: Observer.track_end() return value
    """
    # Initial call to the observer
    if observer is not None:
        observer.track_start(beam)

    # Main loop
    for turn in range(0, turns):
        for i in range(0, line.shape[0]):
            # Symplectic integrators
            if line[i, INDEX_CLASS_CODE] in CLASS_CODE_INTEGRATOR and beam.shape[0]:
                beam = integrators[int(line[i, INDEX_CLASS_CODE])](line[i], beam, **kwargs)
            # Monte-Carlo propagation
            elif line[i, INDEX_CLASS_CODE] in CLASS_CODE_FE and beam.shape[0]:
                beam = mc[int(line[i, INDEX_CLASS_CODE])](line[i], beam, **kwargs)
            # Transfert matrices and tensors
            elif line[i, INDEX_CLASS_CODE] in CLASS_CODE_MATRIX and beam.shape[0]:
                matrix = matrices[int(line[i, INDEX_CLASS_CODE])]
                # For performance considerations, see
                # https://stackoverflow.com/q/48474274/420892
                # Alternative
                # beam = np.einsum('ij,kj->ik', beam, matrix(line[i]))
                beam = beam.dot(matrix(line[i]).T)
            beam = aperture_check(beam, line[i])

            # Per element observation
            if observer is not None and observer.element_by_element_is_active is True:
                if observer.elements is None:  # call the observer for each element
                    observer.element_by_element(turn, i, beam)
                elif i in observer.elements:  # call the observer for given elements
                    observer.element_by_element(turn, i, beam)
        # Per turn observation
        if observer is not None and observer.turn_by_turn_is_active is True:
            observer.turn_by_turn(turn, i, beam)
    # Final call to the observer
    if observer is not None:
        return observer.track_end(turn, i, beam)
    return


def track2(line, beam, turns=1, observer=None, **kwargs):
    """
    Tracking through a beamline.
    Code optimized for performance.
    :param line: beamline description in Manzoni format
    :param beam: initial beam
    :param turns: number of tracking turns
    :param observer: Observer object to witness and record the tracking data
    :param kwargs: optional parameters
    :return: Observer.track_end() return value
    """
    # Initial call to the observer
    if observer is not None:
        observer.track_start(beam)

    # Main loop
    for turn in range(0, turns):
        for i in range(0, line.shape[0]):
            # Symplectic integrators
            if line[i, INDEX_CLASS_CODE] in CLASS_CODE_INTEGRATOR and beam.shape[0]:
                beam = integrators[int(line[i, INDEX_CLASS_CODE])](line[i], beam, **kwargs)
            # Monte-Carlo propagation
            elif line[i, INDEX_CLASS_CODE] in CLASS_CODE_FE and beam.shape[0]:
                beam = mc[int(line[i, INDEX_CLASS_CODE])](line[i], beam, **kwargs)
            # Transfert matrices and tensors
            elif line[i, INDEX_CLASS_CODE] in CLASS_CODE_MATRIX and beam.shape[0]:
                ts = tensors[int(line[i, INDEX_CLASS_CODE])]
                if ts is None:
                    # For performance considerations, see
                    # https://stackoverflow.com/q/48474274/420892
                    # Alternative
                    # beam = np.einsum('ij,kj->ik', beam, matrix(line[i]))
                    r = matrices[int(line[i, INDEX_CLASS_CODE])](line[i], True)  # DO multiply
                    beam = beam.dot(r.T)
                else:
                    # For performance considerations, see
                    # https://stackoverflow.com/q/51046917/420892
                    # Alternatives
                    # np.einsum('ijk,lj,lk->li', T, beam, beam)
                    # ((beam[:, None, None, :]@T).squeeze()@x[..., None]).squeeze()
                    # np.einsum('ijl,lj->li', T@beam.T, beam)
                    rs = matrices[int(line[i, INDEX_CLASS_CODE])](line[i], False)  # DO NOT multiply
                    for u in range(len(rs)):
                        beam = beam.dot(rs[i].T) + np.einsum('lj,ijk,lk->li', beam, ts[i], beam)
            beam = aperture_check(beam, line[i])

            # Per element observation
            if observer is not None and observer.element_by_element_is_active is True:
                if observer.elements is None:  # call the observer for each element
                    observer.element_by_element(turn, i, beam)
                elif i in observer.elements:  # call the observer for given elements
                    observer.element_by_element(turn, i, beam)
        # Per turn observation
        if observer is not None and observer.turn_by_turn_is_active is True:
            observer.turn_by_turn(turn, i, beam)

    # Final call to the observer
    if observer is not None:
        return observer.track_end(turn, i, beam)
    return
