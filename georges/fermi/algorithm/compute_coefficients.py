import os
import sys

import georges_core.kinematics as gkin
import numpy as np
import pandas as pd
import pybdsim
from compute_quantiles import compute_quantile_2dim, compute_quantile_4dim
from lmfit import Model, Parameters

from georges import ureg as _ureg


def cut_data(filepath: str = None, epos: float = 100, matname: str = None, nprimary: int = 10) -> pd.DataFrame:
    """

    Args:
        filepath (str): Path to the file
        epos (float): Theoritical energy at the exit of the degrader
        matname (str): Name of the material
        nprimary (int): Number of particles launched in BDSIM.

    Returns:

    """

    element = "DEG"
    coordinates = ["x", "y", "xp", "yp", "energy", "p", "partID", "parentID", "S"]

    data = pybdsim.Analysis.BDSimOutput(filepath).event.samplers[element].df
    data = data[coordinates]

    data.reset_index(drop=True, inplace=True)

    data.query("parentID == 0 and partID == 2212", inplace=True)  # only get the primaries
    data["energy"] = gkin.etot_to_ekin(data["energy"].values * _ureg.GeV).m_as("MeV")
    data["momentum"] = data["p"] * 1000  # In Mev/c
    data["xp"] *= 1000
    data["yp"] *= 1000
    data["x"] *= 1000
    data["y"] *= 1000

    # Cut the data
    if epos < 220:
        quantile, data_cutted = compute_quantile_4dim(data)

    else:
        quantile, data_cutted = compute_quantile_2dim(data)

    # Properties of the cutted data
    momentum_cut = data_cutted["momentum"].mean()
    deviation_cut = (data_cutted["momentum"] - momentum_cut) / momentum_cut
    dpp_mean_cut = np.mean(deviation_cut)
    dpp_rms_cut = np.std(deviation_cut)
    transmission_cut = len(data_cutted) / nprimary

    # Return the results
    return pd.DataFrame(
        data={
            "MatName": [matname],
            "Epost": [epos],
            "DPPmean_cut": [dpp_mean_cut],
            "DPPrms_cut": [dpp_rms_cut],
            "Transmission_cut": [transmission_cut],
        },
    )


def fit_values(x, **params):
    val = 0.0
    parnames = sorted(params.keys())
    for i, pname in enumerate(parnames):
        val += params[pname] * x**i
    return val


if __name__ == "__main__":
    path = sys.argv[1]
    nparticles = int(sys.argv[2])

    results = pd.DataFrame()
    for file in os.listdir(path):
        if "parquet" in file:
            print(f"Process {file}")
            results = pd.concat(
                [
                    results,
                    cut_data(
                        filepath=os.path.join(path, file),
                        epos=float(file.split("-")[2].replace("E", "").replace("root", "")),
                        matname=file.split("-")[1],
                        nprimary=nparticles,
                    ),
                ],
                axis=0,
            )

    params_transmission = Parameters()
    params_transmission.add("C0", value=1)
    params_transmission.add("C1", value=1)
    params_transmission.add("C2", value=1)
    params_transmission.add("C3", value=0)

    fit_data = Model(fit_values)
    data_transmission = results[["Epost", "Transmission_cut"]]
    fit_transmission = fit_data.fit(
        data_transmission["Transmission_cut"],
        params_transmission,
        x=data_transmission["Epost"],
    )
    coeff_transmission = fit_transmission.best_values

    params_dpp = Parameters()
    params_dpp.add("C0", value=1)
    params_dpp.add("C1", value=1)
    params_dpp.add("C2", value=1)
    params_dpp.add("C3", value=0)
    params_dpp.add("C4", value=0)
    params_dpp.add("C5", value=0)
    data_dpp = results[["Epost", "DPPrms_cut"]]
    fit_dpp = fit_data.fit(data_dpp["DPPrms_cut"], params_dpp, x=data_dpp["Epost"])
    coeff_dpp = fit_dpp.best_values

    print(
        f"""
    materiaName,
    {coeff_transmission['C0']},{coeff_transmission['C1']},{coeff_transmission['C2']},{coeff_transmission['C3']}
    {coeff_dpp['C0']},{coeff_dpp['C1']},{coeff_dpp['C2']},{coeff_dpp['C3']},{coeff_dpp['C4']},{coeff_dpp['C5']}
    """,
    )
