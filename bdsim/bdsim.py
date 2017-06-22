import shutil
import subprocess as sub
import jinja2
import re
import pandas as pd
from .grammar import bdsim_syntax

SUPPORTED_PROPERTIES = ['ANGLE', 'APERTYPE', 'E1', 'E2', 'FINT', 'HGAP', 'THICK', 'TILT']


def element_to_bdsim(e):
    """Convert a pandas.Series representation onto a MAD-X sequence element."""
    bdsim = ""
    if e.KEYWORD in ['MARKER', 'INSTRUMENT']:
        bdsim = "{}: {};".format(e.name.replace('$', ''), "marker")
    if e.KEYWORD in ['DRIFT', 'QUADRUPOLE', 'RBEND', 'SBEND']:
        bdsim = "{}: {}, l={}*m".format(e.name.replace('$', ''), e.KEYWORD.lower(), e.L)
        if pd.notnull(e['ANGLE']):
            bdsim += ", angle=-{}/DEGREE".format(e['ANGLE'])
        #if pd.notnull(e['APERTYPE']):
        #    bdsim += ", aperture={}*m".format(str(e['APERTURE']).strip('[]'))
        bdsim += ';'
    return bdsim


def sequence_to_bdsim(sequence):
    """Convert a pandas.DataFrame sequence onto a BDSim input."""
    sequence.sort_values(by='AT_CENTER', inplace=True)
    if sequence is None:
        return ""
    input = "BRHO=2.3114; DEGREE=pi/180.0;"
    input += " ".join(
        sequence.reset_index().drop_duplicates(subset='index', keep='last').set_index('index').apply(
            element_to_bdsim, axis=1))

    input += "{}: line = ({});".format("ess", ",".join(sequence.index.map(lambda x: x.replace('$', ''))))
    if 'CIRCUIT' in sequence:
        input += '\n'.join(sequence['CIRCUIT'].dropna().map(lambda c: "{}:={{{{ {} or '0.0' }}}};".format(c, c)))
        input += '\n'
    return input


class BDSimException(Exception):
    """Exception raised for errors in the BDSim module."""

    def __init__(self, m):
        self.message = m


class BDSim:
    """A Python wrapper around the BDSim software.

    Sequence and command will be converted with the BDSim grammar and pipe'd to the subprocess.
    """
    def __init__(self, **kwargs):
        self.__input = ""
        self.__beamlines = kwargs.get('beamlines', [])
        self.__path = kwargs.get('path', ".")
        self.__bdsim = kwargs.get('bdsim', None)
        self.__context = kwargs.get('context', {})

        self.__warnings = []
        self.__fatals = []
        self.__output = ""
        self.__template_input = None
        # Convert all sequences to BDSim sequences
        map(self.attach, self.__beamlines)

    def __get_bdsim_path(self):
        return self.__bdsim if self.__bdsim is not None else shutil.which("bdsim")

    def __add_input(self, keyword, strings=()):
        self.__input += madx_syntax[keyword].format(*strings) + '\n'

    def attach(self, beamline):
        self.__beamlines.append(beamline)
        self.__input = sequence_to_bdsim(beamline.line)

    def run(self, context):
        """Run bdsim as a subprocess."""
        self.__template_input = jinja2.Template(self.__input).render(context)
        if self.__get_bdsim_path() is None:
            raise MadxException("Can't run MADX if no valid path and executable are defined.")
        p = sub.Popen([self.__get_bdsim_path()],
                      stdin=sub.PIPE,
                      stdout=sub.PIPE,
                      stderr=sub.STDOUT,
                      cwd=self.__path,
                      shell=True
                      )
        self.__output = p.communicate(input=self.__template_input.encode())[0].decode()
        self.__warnings = [line for line in self.__output.split('\n') if re.search('warning|fatal', line)]
        self.__fatals = [line for line in self.__output.split('\n') if re.search('fatal', line)]
        return self

    def print_input(self):
        """Print the rendered MAD-X input."""
        print(jinja2.Template(self.__input).render(self.context))

    @property
    def path(self):
        """Current MAD-X path."""
        return self.__path

    @property
    def warnings(self):
        """Return warnings from the previous execution run."""
        return self.__warnings

    @property
    def fatals(self):
        """Return fatal errors from the previous execution run."""
        return self.__fatals

    @property
    def input(self):
        """Return the current MAD-X input string."""
        return self.__input

    @property
    def output(self):
        """Return the output of the last MAD-X run."""
        return self.__output

    @property
    def beamlines(self):
        """Return the list of beamlines attached to the instance of MAD-X."""
        return self.__beamlines

    @property
    def context(self):
        """The current state of the beamline."""
        return self.__context

    @context.setter
    def context(self, c):
        self.__context = c

    def print_warnings(self):
        """Print warnings from the previous execution run."""
        [print(w) for w in self.__warnings]

