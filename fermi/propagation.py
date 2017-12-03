from .. import Beamline
from .fermi_eyges import compute_fermi_eyges
from .stopping import residual_energy


class FermiPropagateException(Exception):
    """Exception raised for errors in the Fermi propagate module."""

    def __init__(self, m):
        self.message = m


def propagate(line, beam, db):
    def compute_fermi_eyges_on_slab():
        fe = pd.Series(compute_fermi_eyges(
            db=db,
            material=e['MATERIAL'],
            energy=e['ENERGY'],
            thickness=e['LENGTH'] * 100,
            T=georges.fermi.ICRUProtons
        ))
        fe['NAME'] = e['NAME']
        fe['A0'] = fe['A'][0]
        fe['A1'] = fe['A'][1]
        fe['A2'] = fe['A'][2]
        fe.drop(['A'], inplace=True)
        return fe

    # Do not modify the input beamline
    line_fermi = line.line.copy()

    # Compute energy loss along the line
    energy = beam.energy
    line_fermi.loc[0, "ENERGY"] = energy
    for i, e in line.iterrows():
        if e["TYPE"] == 'slab':
            energy = residual_energy(e['MATERIAL'], e['LENGTH'] * 100, energy, db=db)
            line_fermi.loc[i + 1, "ENERGY"] = energy
        else:
            energy = e["ENERGY"]

    # Beam spreading and scattering following the Fermi-Eyges model
    line_fermi.merge(line_fermi.query("TYPE == 'slab'").apply(compute_fermi_eyges_on_slab, axis=1), how='left')

    # Reconstruct a new beamline
    return Beamline(line_fermi)
