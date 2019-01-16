from typing import Optional, List, Dict
import logging
import numpy as np
from .constants import *
from .. import fermi
from .. import physics
from .. import model as _model
from .. import Beamline
from .. import Beam

FERMI_DB: fermi.MaterialsDB = fermi.MaterialsDB()


def convert_line(line: Beamline,
                 context: Optional[Dict] = None,
                 to_numpy: bool = True,
                 fermi_params: Optional[Dict] = None
                 ):
    context = context or {}
    fermi_params = fermi_params or {}

    def class_conversion(e):
        if e['CLASS'] in ('RFCAVITY',):
            e['CLASS_CODE'] = CLASS_CODES['DRIFT']
        if e['CLASS'] not in CLASS_CODES:
            e['CLASS_CODE'] = CLASS_CODES['NONE']
        else:
            e['CLASS_CODE'] = CLASS_CODES[e['CLASS']]
        return e

    def circuit_conversion(e):
        if 'PLUG' not in e.index.values:
            return e
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
                if context.get(s[0]) is not None:
                    e['APERTURE'] = float(context.get(s[0]))
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
        if e['CLASS'] != 'DEGRADER' and e['CLASS'] != 'SCATTERER':
            return e
        material = e['MATERIAL']
        if str(material) == '' or str(material) == 'vacuum' or e['LENGTH'] == 0:
            return e
        fe = fermi.compute_fermi_eyges(material=str(material),
                                       energy=e['ENERGY_IN'],
                                       thickness=100*e['LENGTH'],
                                       db=FERMI_DB,
                                       t=fermi.DifferentialMoliere,
                                       with_dpp=fermi_params.get('with_dpp', True),
                                       with_losses=fermi_params.get('with_losses', True),
                                       )
        e['FE_A0'] = fe['A'][0]
        e['FE_A1'] = fe['A'][1]
        e['FE_A2'] = fe['A'][2]
        e['FE_DPP'] = fe['DPP']
        e['FE_LOSS'] = fe['LOSS']
        return e

    # Create or copy missing columns
    line_copy = line.copy()
    if 'CLASS' not in line_copy and 'KEYWORD' in line_copy:
        line_copy['CLASS'] = line_copy['KEYWORD']
    if 'CLASS' not in line_copy and 'TYPE' in line_copy:
        line_copy['CLASS'] = line_copy['KEYWORD']

    # Fill with zeros
    for i in INDEX.keys():
        if i not in line_copy:
            line_copy[i] = 0.0
    # Perform the conversion
    line_copy = line_copy.apply(class_conversion, axis=1)
    line_copy = line_copy.apply(circuit_conversion, axis=1)
    line_copy = line_copy.apply(apertype_conversion, axis=1)

    # Energy tracking
    if 'BRHO' not in line.columns:
        if context.get('ENERGY'):
            energy = context['ENERGY']
        elif context.get('PC'):
            energy = physics.momentum_to_energy(context['PC'])
        elif context.get('BRHO'):
            energy = physics.brho_to_energy(context['BRHO'])
        fermi.track_energy(energy, line_copy, FERMI_DB)
        line_copy['BRHO'] = physics.energy_to_brho(line_copy['ENERGY_IN'])
    logging.info(f"Energies: max = {line_copy['ENERGY_IN'].max()}, "
                 f"min = {line_copy['ENERGY_IN'].min()}, "
                 f"loss = {line_copy['ENERGY_IN'].max() - line_copy['ENERGY_IN'].min()}"
                 )

    # Compute Fermi-Eyges parameters
    line_copy = line_copy.apply(fermi_eyges_computations, axis=1)

    # Adjustments for the final format
    line_copy = line_copy.fillna(0.0)
    if to_numpy:
        return line_copy[list(INDEX.keys())].values
    else:
        return line_copy[list(INDEX.keys())+['ENERGY_IN', 'ENERGY_OUT']]


def transform_variables(line, variables) -> list:
    ll = line.reset_index()

    def transform(v):
        i = ll[ll['NAME'] == v[0]].index.values[0]
        j = INDEX[v[1]]
        return [i, j]
    return list(map(transform, variables))


def adjust_line(line, variables, parameters):
    if len(parameters) == 0:
        return line
    it = np.nditer(parameters, flags=['f_index'])
    while not it.finished:
        line[variables[it.index][0], variables[it.index][1]] = it[0]
        it.iternext()
    return line


def transform_elements(line, elements: List[str]) -> List[int]:
    """Convert a list of elements (referenced by names) onto a list of indices"""
    ll = line.reset_index()

    def transform(e):
        return ll[ll['NAME'] == e].index.values[0]
    return list(map(transform, elements))


def _process_model_argument(model,
                            line: Optional[Beamline],
                            beam: Optional[Beam],
                            context: Optional[Dict],
                            exception: type = Exception,
                            ) -> dict:
    """

    :param model:
    :param line:
    :param beam:
    :param context:
    :param exception:
    :return:
    """
    if model is None:
        if line is None or beam is None:
            raise exception("Beamline and Beam objects need to be defined.")
        else:
            # Beamline
            if isinstance(line, Beamline):
                georges_line = line
                manzoni_line = convert_line(line.line, context)
            else:
                raise exception("'line' must be a Georges Beamline")
            # Beam
            if isinstance(beam, Beam):
                georges_beam = beam
                manzoni_beam = np.array(beam.distribution)
            else:
                raise exception("'line' must be a Georges Beam")
            # Context
            georges_context = context
    else:
        # Access the object's model
        if (not isinstance(model, _model.Model) and not isinstance(model, _model.ManzoniModel)) \
                and hasattr(model, 'model'):
            model = model.model

        if isinstance(model, _model.ManzoniModel):
            georges_line = model.gbeamline
            georges_beam = model.gbeam
            georges_context = model.context
            manzoni_line = model.beamline
            manzoni_beam = model.beam
        elif isinstance(model, _model.Model):
            georges_line = model.beamline
            georges_beam = model.beam
            georges_context = model.context
            manzoni_model = _model.ManzoniModel(model)
            manzoni_line = manzoni_model.beamline
            manzoni_beam = manzoni_model.beam
        else:
            raise exception("'model' must be a Georges Model or a Georges ManzoniModel.")

    return {
        'georges_context': georges_context,
        'georges_line': georges_line,
        'georges_beam': georges_beam,
        'manzoni_line': manzoni_line,
        'manzoni_beam': manzoni_beam,
    }
