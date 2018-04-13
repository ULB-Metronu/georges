import numpy as np
from .transfer import transfer
from .kick import kick
from .constants import *
from .aperture import aperture_check
from .. import fermi
from .. import physics


def convert_line(line, context={}, to_numpy=True):
    def class_conversion(e):
        if e['CLASS'] in ('RFCAVITY', 'HKICKER'):
            e['CLASS_CODE'] = CLASS_CODES['DRIFT']
        if e['CLASS'] not in CLASS_CODES:
            e['CLASS_CODE'] = CLASS_CODES['NONE']
        else:
            e['CLASS_CODE'] = CLASS_CODES[e['CLASS']]
        return e

    def circuit_conversion(e):
        if e['PLUG'] in INDEX and e['PLUG'] != 'APERTURE':
            e[e['PLUG']] = context.get(e['CIRCUIT'], 0.0)
        return e

    def apertype_conversion(e):
        # Default aperture
        if 'APERTYPE' not in e.index.values:
            e['APERTYPE_CODE'] = APERTYPE_CODE_NONE
            e['APERTURE'] = 0.0
            e['APERTURE_2'] = 0.0
            return e
        # Aperture types
        if e['APERTYPE'] == 'CIRCLE':
            e['APERTYPE_CODE'] = APERTYPE_CODE_CIRCLE
        elif e['APERTYPE'] == 'RECTANGLE':
            e['APERTYPE_CODE'] = APERTYPE_CODE_RECTANGLE
        else:
            e['APERTYPE_CODE'] = APERTYPE_CODE_NONE
            e['APERTURE'] = 0.0
            e['APERTURE_2'] = 0.0
        # Aperture sizes
        if not isinstance(e['APERTURE'], str):
            if np.isnan(e['APERTURE']) and e['PLUG'] == 'APERTURE':
                s = e['CIRCUIT'].strip('[{}]').split(',')
                if context.get(s[0], 0):
                    e['APERTURE'] = float(context.get(s[0], 1.0))
                else:
                    e['APERTURE'] = 1.0
                if len(s) > 1:
                    e['APERTURE_2'] = float(context.get(s[1], 1.0))
                else:
                    e['APERTURE_2'] = 1.0
        else:
            s = e['APERTURE'].strip('[{}]').split(',')
            e['APERTURE'] = float(s[0])
            if len(s) > 1:
                e['APERTURE_2'] = float(s[1])
        return e

    def fermi_eyges_computations(e):
        if e['CLASS'] != 'DEGRADER':
            return e
        fe = fermi.compute_fermi_eyges(material=e['MATERIAL'],
                                       energy=e['ENERGY_IN'],
                                       thickness=100*e['LENGTH'],
                                       db=db,
                                       t=fermi.DifferentialMoliere)
        e['FE_A0'] = fe['A'][0]
        e['FE_A1'] = fe['A'][1]
        e['FE_A2'] = fe['A'][2]
        e['FE_DPP'] = 0
        return e

    # Create or copy missing columns
    line = line.copy()
    if 'CLASS' not in line and 'KEYWORD' in line:
        line['CLASS'] = line['KEYWORD']
    if 'CLASS' not in line and 'TYPE' in line:
        line['CLASS'] = line['KEYWORD']

    # Fill with zeros
    for i in INDEX.keys():
        if i not in line:
            line[i] = 0.0

    # Perform the conversion
    line = line.apply(class_conversion, axis=1)
    line = line.apply(circuit_conversion, axis=1)
    line = line.apply(apertype_conversion, axis=1)

    # Energy tracking
    db = fermi.MaterialsDB()
    energy = 230  # Default value
    if context.get('ENERGY'):
        energy = context['ENERGY']
    elif context.get('PC'):
        energy = physics.momentum_to_energy(context['PC'])
    fermi.track_energy(energy, line, db)
    line['BRHO'] = physics.energy_to_brho(line['ENERGY_IN'])

    # Compute Fermi-Eyges parameters
    line = line.apply(fermi_eyges_computations, axis=1)

    # Adjustments for the final format
    line = line.fillna(0.0)
    if to_numpy:
        return line[list(INDEX.keys())].as_matrix()
    else:
        print(CLASS_CODES)
        return line[list(INDEX.keys())]


def transform_variables(line, variables):
    ll = line.reset_index()

    def transform(v):
        i = ll[ll['index'] == v[0]].index.values[0]
        j = INDEX[v[1]]
        return [i, j]
    return list(map(transform, variables))


def adjust_line(line, variables, parameters):
    for i, p in enumerate(parameters):
        line[variables[i][0], variables[i][1]] = p
    return line  # Check later if calling this function will create a copy


def track(line, beam, turns=1, observer=None, **kwargs):
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
    for j in range(0, turns):
        for i in range(0, line.shape[0]):
            if line[i, INDEX_CLASS_CODE] in CLASS_CODE_KICK and beam.shape[0] > 0:
                # In place operation
                beam = kick[int(line[i, INDEX_CLASS_CODE])](line[i], beam, **kwargs)
            elif line[i, INDEX_CLASS_CODE] in CLASS_CODE_MATRIX:
                matrices = transfer[int(line[i, INDEX_CLASS_CODE])]
                # For performance considerations, see
                # https://stackoverflow.com/q/48474274/420892
                # b = np.einsum('ij,kj->ik', b, matrix(line[i]))
                beam = beam.dot(matrices(line[i]).T)
            beam = aperture_check(beam, line[i])
            if observer is not None and observer.element_by_element_is_active is True:
                observer.element_by_element(j, i, beam)
        if observer is not None and observer.turn_by_turn_is_active is True:
            observer.turn_by_turn(j, i, beam)
    # Final call to the observer
    if observer is not None:
        return observer.track_end(j, i, beam)
    return
