import pandas as pd
import numpy as np
from .. import Beamline
from .fermi_eyges import compute_fermi_eyges
from .stopping import residual_energy
from .mcs import DifferentialMoliere
from .materials import Vacuum


class FermiPropagateException(Exception):
    """Exception raised for errors in the Fermi propagate module."""

    def __init__(self, m):
        self.message = m


def track_energy(energy, line_fermi, db):
    for i, e in line_fermi.iterrows():
        line_fermi.loc[i, 'ENERGY_IN'] = energy
        if e['TYPE'] == 'slab' or e['TYPE'] == 'gap' or e['CLASS'] == 'DEGRADER':
            energy = residual_energy(
                e['MATERIAL'], e['LENGTH'] * 100, energy, db=db
            ) if str(e['MATERIAL']) != 'vacuum' else energy
            line_fermi.loc[i, 'ENERGY_OUT'] = energy
            line_fermi.loc[i, 'DeltaE'] = line_fermi.loc[i, 'ENERGY_IN'] - energy
        else:
            line_fermi.loc[i, 'ENERGY_OUT'] = energy


def propagate(line, beam, db, model=DifferentialMoliere, gaps='vacuum'):
    def compute_fermi_eyges_on_slab(slab):
        # If vacuum, return a null-element
        if str(slab['MATERIAL']) == Vacuum:
            return pd.Series({
                'A0': 0,
                'A1': 0,
                'A2': 0,
                'B': 0
            }).rename(e.name)
        # Do the actual computation with other materials.py
        fe = compute_fermi_eyges(
            db=db,
            material=str(slab['MATERIAL']),
            energy=slab['ENERGY_IN'],
            thickness=slab['LENGTH'] * 100,
            t=model
        )
        return pd.Series({
            'A0': fe['A'][0],
            'A1': fe['A'][1],
            'A2': fe['A'][2],
            'B': fe['B']
        }).rename(slab.name)

    # Default beam
    if beam is None:
        beam = {
            'A0': 0,
            'A1': 0,
            'A2': 0,
        }

    # Do not modify the input beamline
    line_fermi = line.line.copy()

    # Add air gaps
    if gaps is not None:
        air_gaps = pd.DataFrame()
        air_gaps['AT_ENTRY'] = line_fermi['AT_EXIT'].shift(1)
        air_gaps['MATERIAL'] = gaps
        air_gaps['LENGTH'] = (line_fermi['AT_ENTRY'] - line_fermi['AT_EXIT'].shift(1)).values
        air_gaps['NAME'] = np.roll(line_fermi.index.values, 1) + '_' + air_gaps.index + '_gap'
        air_gaps['TYPE'] = 'gap'
        air_gaps.dropna(inplace=True, subset=['LENGTH'])
        air_gaps.set_index('NAME', inplace=True)
        with_gaps = pd.concat([line_fermi, air_gaps], sort=False).sort_values(by='AT_ENTRY')
        line_fermi = with_gaps

    # Add columns as needed
    if 'TYPE' not in line_fermi:
        line_fermi['TYPE'] = line_fermi['CLASS']

    # Compute energy loss along the line
    track_energy(beam['energy'], line_fermi, db)

    # Beam spreading and scattering following the Fermi-Eyges model
    line_fermi = line_fermi.merge(
        line_fermi.query("TYPE == 'slab' or TYPE == 'gap'").apply(compute_fermi_eyges_on_slab, axis=1),
        how='left',
        left_index=True,
        right_index=True,
        )

    # Sum the matrix elements
    acc = [beam['A0'], beam['A1'], beam['A2']]
    line_fermi['LENGTH'].fillna(0.0, inplace=True)
    line_fermi['A0'].fillna(0.0, inplace=True)
    line_fermi['A1'].fillna(0.0, inplace=True)
    line_fermi['A2'].fillna(0.0, inplace=True)
    line_fermi['B'].fillna(0.0, inplace=True)
    for i, e in line_fermi.iterrows():
        line_fermi.loc[i, 'A0_IN'] = acc[0]
        line_fermi.loc[i, 'A1_IN'] = acc[1]
        line_fermi.loc[i, 'A2_IN'] = acc[2]
        line_fermi.loc[i, 'A0_OUT'] = acc[0] + e['A0']
        line_fermi.loc[i, 'A1_OUT'] = acc[0] * e['LENGTH'] + acc[1] + e['A1']
        line_fermi.loc[i, 'A2_OUT'] = acc[0] * e['LENGTH'] ** 2 + 2 * acc[1] * e['LENGTH'] + acc[2] + e['A2']
        acc = [line_fermi.loc[i, 'A0_OUT'], line_fermi.loc[i, 'A1_OUT'], line_fermi.loc[i, 'A2_OUT']]
    line_fermi['B_IN'] = line_fermi['A0_IN'] * line_fermi['A2_IN'] - line_fermi['A1_IN'] ** 2
    line_fermi['B_OUT'] = line_fermi['A0_OUT'] * line_fermi['A2_OUT'] - line_fermi['A1_OUT'] ** 2

    return Beamline(line_fermi)
