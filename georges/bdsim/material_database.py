

def get_bdsim_material(mat_name):

    if mat_name == 'tantalum':
        return 'G4_Ta'

    if mat_name == 'lexan':
        return 'G4_POLYCARBONATE'

    if mat_name == 'polyethylen':
        return 'G4_POLYETHYLENE'

    if mat_name == 'brass':
        return 'G4_BRASS'
