"""
TODO
"""
from typing import Mapping, Optional

import pandas as _pd
from georges_core.sequences import Sequence

from .. import ureg as _ureg
from .mcs import DifferentialMoliere as _DifferentialMoliere


def track_energy(sequence: Sequence, energy: _ureg.Quantity):
    s: _pd.DataFrame = sequence.to_df()
    energy = energy.m_as("MeV")
    if "MATERIAL" not in s.columns:
        s["ENERGY_IN"] = energy
        s["ENERGY_OUT"] = energy
        return s
    for i, e in s.iterrows():
        s.loc[i, "ENERGY_IN"] = energy
        if not _pd.isnull(e["MATERIAL"]) and e["L"] != 0.0:
            e = e["MATERIAL"].stopping(thickness=e["L"], kinetic_energy=energy * _ureg.MeV)
            if e is not None:
                energy = e.ekin.m_as("MeV")
        s.loc[i, "ENERGY_OUT"] = energy
        s.loc[i, "DELTA_E"] = s.at[i, "ENERGY_IN"] - energy
    return s


def propagate(sequence: Sequence, energy: _ureg.Quantity, beam: Optional[Mapping] = None, model=_DifferentialMoliere):
    """

    Args:
        sequence: Sequence to track in
        energy: Beam energy
        beam: A_i beam parameters
        model: Model to use. Default is DifferentialMoliere

    Returns:

    """
    # Default beam
    if beam is None:
        beam = {
            "A0": 0,
            "A1": 0,
            "A2": 0,
        }

    # Compute energy loss along the line
    s = track_energy(sequence=sequence, energy=energy)

    # Beam spreading and scattering following the Fermi-Eyges model
    for i, e in s.iterrows():
        if not _pd.isnull(e["MATERIAL"]):
            fe = e["MATERIAL"].scattering(
                kinetic_energy=energy,
                thickness=e["L"],
                model=model,
            )
        else:
            fe = {"A": (0, 0, 0), "B": 0}
        s.at[i, "A0"] = fe["A"][0]
        s.at[i, "A1"] = fe["A"][1]
        s.at[i, "A2"] = fe["A"][2]
        s.at[i, "B"] = fe["B"]

    # Sum the matrix elements
    acc = [beam["A0"], beam["A1"], beam["A2"]]
    for i, e in s.iterrows():
        s.loc[i, "A0_IN"] = acc[0]
        s.loc[i, "A1_IN"] = acc[1]
        s.loc[i, "A2_IN"] = acc[2]
        s.loc[i, "A0_OUT"] = acc[0] + e["A0"]
        s.loc[i, "A1_OUT"] = acc[0] * e["L"].m_as("m") + acc[1] + e["A1"]
        s.loc[i, "A2_OUT"] = acc[0] * e["L"].m_as("m") ** 2 + 2 * acc[1] * e["L"].m_as("m") + acc[2] + e["A2"]
        acc = [s.loc[i, "A0_OUT"], s.loc[i, "A1_OUT"], s.loc[i, "A2_OUT"]]
    s["B_IN"] = s["A0_IN"] * s["A2_IN"] - s["A1_IN"] ** 2
    s["B_OUT"] = s["A0_OUT"] * s["A2_OUT"] - s["A1_OUT"] ** 2

    return s
