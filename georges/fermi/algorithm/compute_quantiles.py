from typing import Any, List, Tuple, Union

import numpy as np
import pandas as pd
from scipy.optimize import minimize


def compute_matrix(beam: pd.DataFrame = None) -> List[Union[np.ndarray, Any]]:
    """Return the transformation matrix"""

    def get_matrix(u, v):
        cov = np.cov(u, v)
        emit = np.sqrt(np.linalg.det(cov))
        beta = cov[0, 0] / emit
        alpha = -cov[0, 1] / emit

        return [[1, 0], [alpha, beta]]

    x = beam["x"].values
    xp = beam["xp"].values
    y = beam["y"].values
    yp = beam["yp"].values

    return [get_matrix(x, xp), get_matrix(y, yp)]


def normalize_phase_space(beam: pd.DataFrame = None, matrix_h: np.ndarray = None, matrix_v: np.ndarray = None):
    bnorm = beam.copy(deep=True)

    b = [beam["x"].values, beam["xp"].values]
    d = [beam["y"].values, beam["yp"].values]

    # Normalize the coordinates
    bnorm["x"], bnorm["xp"] = np.dot(matrix_h, b)
    bnorm["y"], bnorm["yp"] = np.dot(matrix_v, d)

    return bnorm


# 2D algorithm
def compute_quantile_2dim(data: pd.DataFrame = None) -> Tuple[float, pd.DataFrame]:
    """
    Compute the quantile in two dimensions: momentum
    Args:
        data: initial beam dataframe

    Returns:
        optimize_quant: optimised quantile
        beam: optimised beam

    """

    beam = data.copy()

    beam["ptrans"] = beam["xp"] ** 2 + beam["yp"] ** 2

    # Compute the quantile
    quantile = minimize(
        lambda xi: get_deviation_2dim(
            data=beam,
            quant=xi,
        ),
        0,
        method="Nelder-Mead",
        tol=1e-8,
    )

    optimize_quant = quantile.x[0]

    # Apply the quantile
    beam = beam[beam["ptrans"] <= beam["ptrans"].quantile(1 - optimize_quant)]

    return optimize_quant, beam


def compute_chisquare_2dim(
    data_x: np.array = None,
    data_y: np.array = None,
) -> float:
    """Computes the chi-square distribution"""

    u_c = (data_x - data_x.mean()) / data_x.std()
    v_c = (data_y - data_y.mean()) / data_y.std()
    p_c = u_c**2 + v_c**2

    return p_c.std()


def get_deviation_2dim(data: pd.DataFrame = None, quant: np.array = None) -> float:
    data_quantile = data[data["ptrans"] <= data["ptrans"].quantile(1 - quant[0])]
    data_x = data_quantile["xp"].values
    data_y = data_quantile["yp"].values
    deviation_chisquare = np.abs((compute_chisquare_2dim(data_x, data_y)) - 2)
    return deviation_chisquare


def get_cutted_data_2dim(data: pd.DataFrame = None, quantile: float = 0) -> pd.DataFrame:
    beam = data.copy()

    beam["ptrans"] = beam["xp"] ** 2 + beam["yp"] ** 2
    beam = beam[beam["trans"] <= beam["trans"].quantile(1 - quantile)]

    return beam


# 4D algorithm
def compute_quantile_4dim(beam: pd.DataFrame = None) -> Tuple[float, pd.DataFrame]:
    """
    Compute the quantile in 4 dimension
    Args:
        beam: initial beam dataframe

    Returns:
        optimize_quant: optimized quantile
        beam: optimized beam
    """

    # Compute the normalisation matrix
    matrix_h, matrix_v = compute_matrix(beam)

    # Normalize the coordinates
    beam = normalize_phase_space(beam, matrix_h, matrix_v)

    # Compute the transversal coordinates
    beam["trans"] = beam["x"] ** 2 + beam["y"] ** 2 + beam["xp"] ** 2 + beam["yp"] ** 2

    # Compute the quantile
    quantile = minimize(
        lambda xi: get_deviation_4dim(
            data=beam,
            quant=xi,
        ),
        0,
        method="Nelder-Mead",
        tol=1e-8,
        options={"disp": False, "fatol": 1e-12},
    )

    optimize_quant = quantile.x[0]

    # Apply the quantile
    beam = beam[beam["trans"] <= beam["trans"].quantile(1 - optimize_quant)]

    # Denormalize the coordinates, use the inverse of the matrix
    beam = normalize_phase_space(beam, np.linalg.inv(matrix_h), np.linalg.inv(matrix_v))

    return optimize_quant, beam


def compute_chisquare_4dim(
    data_x: np.array = None,
    data_xp: np.array = None,
    data_y: np.array = None,
    data_yp: np.array = None,
) -> float:
    """Computes the chi-square distribution"""

    xn = (data_x - data_x.mean()) / data_x.std()
    xpn = (data_xp - data_xp.mean()) / data_xp.std()

    yn = (data_y - data_y.mean()) / data_y.std()
    ypn = (data_yp - data_yp.mean()) / data_yp.std()

    p_c = xn**2 + xpn**2 + yn**2 + ypn**2
    return p_c.var()


def get_deviation_4dim(data: pd.DataFrame = None, quant=None) -> float:
    if quant is None:
        quant = [0]
    data_quantile = data[data["trans"] <= data["trans"].quantile(1 - quant[0])]

    # Compute the chi-square std
    data_x = data_quantile["x"].values
    data_xp = data_quantile["xp"].values
    data_y = data_quantile["y"].values
    data_yp = data_quantile["yp"].values
    deviation_chisquare = np.abs(
        (
            compute_chisquare_4dim(
                data_x,
                data_xp,
                data_y,
                data_yp,
            )
        )
        - 8.0,
    )
    return deviation_chisquare


def get_cutted_data_4dim(beam: pd.DataFrame = None, quantile: float = 0) -> pd.DataFrame:

    # Compute the normalisation matrix
    matrix_h, matrix_v = compute_matrix(beam)

    # Normalize the coordinates
    beam = normalize_phase_space(beam, matrix_h, matrix_v)

    # Compute the transversal coordinates and apply quantile
    beam["trans"] = beam["x"] ** 2 + beam["y"] ** 2 + beam["xp"] ** 2 + beam["yp"] ** 2
    beam = beam[beam["trans"] <= beam["trans"].quantile(1 - quantile)]

    # Denormalize the coordinates
    beam = normalize_phase_space(beam, np.linalg.inv(matrix_h), np.linalg.inv(matrix_v))

    return beam
