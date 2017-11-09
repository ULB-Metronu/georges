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
    """TODO
    """
    # Process arguments
    line = kwargs.get('line', None)
    if line is None:
        raise MatchException("Beamline and MAD-X objects need to be defined.")
    m = Madx(beamlines=[line])
    m.beam(line.name)
    m.match(sequence=line.name, line=True, vary=['IQ1E', 'IQ2E'],
            constraints={'B1B1': ('BETX', 1.0), 'B1B2': ('ALFX', 1.0)}, context=context)
    errors = m.run(**kwargs).fatals
    if kwargs.get("debug", False):
        print(m.input)
    if len(errors) > 0:
        print(errors)
        raise MatchException("MAD-X ended with fatal error during matching.")
    return process_match_output(m.output)


def process_match_output(output):
    regex_target = re.compile("Final Penalty Function =   (.*)")
    regex_variable = re.compile("(.+)\s+(-?\d\.\d*e.\d\d)\s+(-?\d\.\d*e.\d\d)\s+(-?\d\.\d*e.\d\d)\s+(-?\d\.\d*e.\d\d)")
    regex_constraints = re.compile(
        "(.+):?\d?\s+(betx)\s+(\d)\s+(-?\d\.\d*E.\d\d)\s+(-?\d\.\d*E.\d\d)\s+(-?\d\.\d*E.\d\d)")

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
        'penalty': match_penalty,
        'summary': filter(None, match_summary)
    }
