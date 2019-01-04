import re
import numpy as np
from .madx import Madx

MADX_MATCH_START_SUMMARY_STRING = 'MATCH SUMMARY'
MADX_MATCH_END_SUMMARY_STRING = 'END MATCH SUMMARY'


class MatchException(Exception):
    """Exception raised for errors in the Match module."""

    def __init__(self, m):
        self.message = m


def match(**kwargs):
    """Interface to the MAD-X matchnig module.
    """
    # Process arguments
    line = kwargs.get('line', None)
    if line is None:
        raise MatchException("Beamline and MAD-X objects need to be defined.")
    m = Madx(beamlines=[line], ptc_use_knl_only=kwargs.get('ptc_use_knl_only', False))
    m.beam(line.name)
    m.match(sequence=line.name,
            line=not kwargs.get('periodic', True),
            vary=kwargs.get('vary', {}),
            constraints=kwargs.get('constraints', []),
            global_constraints=kwargs.get('global_constraints', []),
            context=kwargs.get('context', {}),
            method=kwargs.get('method', 'jacobian'),
            ptc=kwargs.get('ptc', False),
            ptc_params=kwargs.get('ptc_params', {})
            )
    errors = m.run(**kwargs).fatals
    if kwargs.get("debug", False):
        print(m.input)
    if len(errors) > 0:
        print(errors)
        raise MatchException("MAD-X ended with fatal error during matching.")
    data = process_match_output(m.output)
    matched_context = kwargs.get('context', {}).copy()
    for k, v in data['variables'].items():
        matched_context[k.upper()] = float(v['final'])
    data['context'] = matched_context
    return data


def process_match_output(output):
    regex_target = re.compile("Final Penalty Function =   (.*)")
    regex_variable = re.compile(r"(.+)\s+(-?\d\.\d*e.\d\d)\s+(-?\d\.\d*e.\d\d)\s+(-?\d\.\d*e.\d\d)\s+(-?\d\.\d*e.\d\d)")
    regex_constraints = re.compile(
        r"(.+):?\d?\s+(betx)\s+(\d)\s+(-?\d\.\d*E.\d\d)\s+(-?\d\.\d*E.\d\d)\s+(-?\d\.\d*E.\d\d)")

    match_flag = False
    match_variables = {}
    match_constraints = {}
    match_penalty = np.inf
    match_summary = []
    for line in output.splitlines():
        if line.startswith(MADX_MATCH_START_SUMMARY_STRING):
            match_flag = True
        if match_flag:
            match_summary.append(line)
            for m in re.finditer(regex_constraints, line):
                v = list(map(str.strip, m.groups()))
                if not match_constraints.get(v[0]):
                    match_constraints[v[0]] = {}
                match_constraints[v[0]][v[1]] = {'target': v[3], 'final': v[4]}
            for m in re.finditer(regex_variable, line):
                v = list(map(str.strip, m.groups()))
                match_variables[v[0]] = {'initial': v[2], 'final': v[1]}
            for m in re.finditer(regex_target, line):
                v = list(map(str.strip, m.groups()))
                match_penalty = v[0]
        if line.startswith(MADX_MATCH_END_SUMMARY_STRING):
            match_flag = False
    return {
        'constraints': match_constraints,
        'variables': match_variables,
        'penalty': float(match_penalty),
        'summary': list(filter(None, match_summary))
    }
