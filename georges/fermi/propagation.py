import pandas as _pd
from . import materials
from .. import ureg as _ureg
from georges_core.sequences import Sequence


def track_energy(sequence: Sequence, energy: _ureg.Quantity):
    s: _pd.DataFrame = sequence.to_df()
    energy = energy.m_as('MeV')
    if 'MATERIAL' not in s.columns:
        s['ENERGY_IN'] = energy
        s['ENERGY_OUT'] = energy
        return s
    for i, e in s.iterrows():
        s.loc[i, 'ENERGY_IN'] = energy
        if not _pd.isnull(e['MATERIAL']) and e['L'] != 0.0:
            m = getattr(materials, e['MATERIAL'].title())
            energy = m.stopping(thickness=e['L'], kinetic_energy=energy * _ureg.MeV).ekin.m_as('MeV')
        s.loc[i, 'ENERGY_OUT'] = energy
        s.loc[i, 'DELTA_E'] = s.at[i, 'ENERGY_IN'] - energy
    return s
