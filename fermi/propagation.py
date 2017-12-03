import pandas as pd
from .. import Beamline
from .fermi_eyges import compute_fermi_eyges
from .stopping import residual_energy
from .mcs import DifferentialMoliere


class FermiPropagateException(Exception):
    """Exception raised for errors in the Fermi propagate module."""

    def __init__(self, m):
        self.message = m


def propagate(line, parameters, db):
    def compute_fermi_eyges_on_slab(slab):
        fe = compute_fermi_eyges(
            db=db,
            material=slab['MATERIAL'],
            energy=slab['ENERGY'],
            thickness=slab['LENGTH'] * 100,
            T=DifferentialMoliere
        )
        return pd.Series({
            'A0': fe['A'][0],
            'A1': fe['A'][1],
            'A2': fe['A'][2]
        }).rename(e.name)

    # Do not modify the input beamline
    line_fermi = line.line.copy()

    # Compute energy loss along the line
    energy = parameters['energy']
    for i, e in line_fermi.iterrows():
        line_fermi.loc[i, 'ENERGY'] = energy
        if e["TYPE"] == 'slab':
            energy = residual_energy(e['MATERIAL'], e['LENGTH'] * 100, energy, db=db)

    # Beam spreading and scattering following the Fermi-Eyges model
    return Beamline(
        line_fermi.merge(line_fermi.query("TYPE == 'slab'").apply(compute_fermi_eyges_on_slab, axis=1),
                         how='left',
                         left_index=True,
                         right_index=True,
                         )
    )
