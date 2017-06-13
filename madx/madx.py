import shutil
import subprocess as sub
import jinja2
import re
import pandas as pd
from .grammar import madx_syntax

SUPPORTED_PROPERTIES = ['ANGLE', 'APERTYPE', 'E1', 'E2', 'FINT', 'HGAP', 'THICK', 'TILT']


def element_to_mad(e):
    """Convert a pandas.Series representation onto a MAD-X sequence element."""
    if not e['PHYSICAL'] or pd.isnull(e['PHYSICAL']):
        return ""
    mad = "{}: {}, ".format(e.name, e.CLASS)
    mad += ', '.join(["{}={}".format(p, e[p]) for p in SUPPORTED_PROPERTIES if pd.notnull(e[p])])
    if pd.notnull(e['ORBIT_LENGTH']):
        mad += ", L={}".format(e['ORBIT_LENGTH'])
    if pd.notnull(e['APERTYPE']):
        mad += ", APERTURE={}".format(str(e['APERTURE']).strip('[]'))
    if pd.notnull(e.get('PLUG')) and pd.notnull(e.get('CIRCUIT')):
        mad += ", {}:={}".format(e['PLUG'], e['CIRCUIT'])
    mad += ", AT={}".format(e['AT_CENTER'])
    mad += ";"
    return mad


def sequence_to_mad(sequence):
    """Convert a pandas.DataFrame sequence onto a MAD-X input."""
    sequence.sort_values(by='AT_CENTER', inplace=True)
    if sequence is None:
        return ""
    input = "{}: SEQUENCE, L={}, REFER=CENTER;\n".format(sequence.name, sequence.length)
    input += '\n'.join(sequence.apply(element_to_mad, axis=1)) + '\n'
    input += "ENDSEQUENCE;\n"
    if 'CIRCUIT' in sequence:
        input += '\n'.join(sequence['CIRCUIT'].dropna().map(lambda c: "{}:={{{{ {} or '0.0' }}}};".format(c, c)))
        input += '\n'
    return input


class MadxException(Exception):
    """Exception raised for errors in the Madx module."""

    def __init__(self, m):
        self.message = m


class Madx:
    """A Python wrapper around the MAD-X executable.

    Sequence and command will be converted with the MAD-X grammar and pipe'd to the subprocess.
    """
    def __init__(self, **kwargs):
        self.__input = ""
        self.__beamlines = kwargs.get('beamlines', [])
        self.__path = kwargs.get('path', ".")
        self.__madx = kwargs.get('madx', None)
        self.__context = kwargs.get('context', {})

        self.__warnings = []
        self.__fatals = []
        self.__output = ""
        self.__template_input = None
        # Convert all sequences to MAD-X sequences
        map(self.attach, self.__beamlines)

    def __get_madx_path(self):
        return self.__madx if self.__madx is not None else shutil.which("madx")

    def __add_input(self, keyword, strings=()):
        self.__input += madx_syntax[keyword].format(*strings) + '\n'

    def attach(self, beamline):
        self.__beamlines.append(beamline)
        self.__input = sequence_to_mad(beamline.line)

    def run(self, context):
        """Run madx as a subprocess."""
        self.__input += madx_syntax['stop']
        self.__template_input = jinja2.Template(self.__input).render(context)
        if self.__get_madx_path() is None:
            raise MadxException("Can't run MADX if no valid path and executable are defined.")
        p = sub.Popen([self.__get_madx_path()],
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

    def set(self, k, v):
        """Set a single variable in the context. Allows method chaining."""
        self.__context[k] = v
        return sel

    def print_warnings(self):
        """Print warnings from the previous execution run."""
        [print(w) for w in self.__warnings]

    def raw(self, raw):
        """Add a raw MAD-X command to the input."""
        self.__input += raw + "\n"
        return self

    def select_columns(self, flag, columns):
        """Add a MAD-X `select` command."""
        self.__add_input('select_columns', (flag, *columns))
        return self

    def call_file(self, file):
        """Add a MAD-X `call` command."""
        self.__add_input('call_file', (file,))
        return self

    def beam(self, line_name):
        """Add a MAD-X `beam` command."""
        self.__add_input('beam')
        self.use_sequence(line_name)
        return self

    def use_sequence(self, sequence):
        """Add a MAD-X `use sequence` command."""
        self.__add_input('use_sequence', (sequence,))
        return self

    def rbarc(self):
        """Add a (legacy) MAD-X `rbarc` option."""
        self.__add_input('rbarc')
        return self

    def twiss(self, **kwargs):
        """Add a (ptc) `twiss` MAD-X command."""
        if kwargs.get('ptc'):
            self.__ptc_twiss(**kwargs)
        else:
            self.__twiss(**kwargs)

    def __ptc_twiss(self, **kwargs):
        self.__add_input('ptc_create_universe')
        self.__add_input('ptc_create_layout',
                         (False, 1, 4, 4, True))
        self.__add_input('ptc_twiss_beamline', (kwargs.get('file', 'ptc_twiss.outx'),))
        self.__add_input('ptc_end')

    def __twiss(self, **kwargs):
        options = ""
        for k, v in kwargs.items():
            if k not in ['ptc']:
                options += ",%s=%s" % (k,v)
        self.__add_input('twiss_beamline', (kwargs.get('file', 'twiss.outx'), options))
        return self

    def makethin(self, sequence, **kwargs):
        """Add a MAD-X `makethin` command."""
        style = kwargs.get('style', 'TEAPOT')
        dipole_slices = kwargs.get('dipole_slices', 4)
        quadrupole_slices = kwargs.get('quadrupole_slices', 4)
        self.__input += "SELECT, FLAG=makethin, CLASS=quadrupole, THICK=false, SLICE={};\n".format(quadrupole_slices)
        self.__input += "SELECT, FLAG=makethin, CLASS=rbend, THICK=false, SLICE={};\n".format(dipole_slices)
        self.__add_input('makethin', (sequence, style))
        self.use_sequence(sequence)
        return self

    def __add_particles_for_tracking(self, particles, ptc=False):
        if {'X', 'PX', 'Y', 'PY', 'DPP'} > set(particles):
            return
        for r in particles.iterrows():
            if ptc:
                self.__add_input('ptc_start', tuple(r[1]))
            else:
                self.__add_input('start_particle', tuple(r[1]))

    def match(self, **kwargs):
        seq = kwargs.get("sequence", None)
        vary = kwargs.get("vary", None)
        constraints = kwargs.get("constraints", None)
        if seq is None:
            raise MadxException("A sequence name must be provided.")
        if vary is None or len(vary) < 1:
            raise MadxException("A list of length > 0 of parameters must be provided.")
        if constraints is None:
            raise MadxException("A dictionary of constraints should be provided.")
        self.__add_input('match', (sequence,))

    def track(self, particles, beamline, **kwargs):
        """Add a ptc `track` command."""
        if kwargs.get('ptc', True):
            self.__ptc_track(particles, beamline, **kwargs)
        else:
            self.__track(particles, beamline, **kwargs)

    def __track(self, particles, beamline, **kwargs):
        if len(particles) == 0:
            print("No particles to track... Doing nothing.")
            return
        self.makethin(beamline.name, **kwargs)
        self.__add_input('track_beamline')
        self.__add_particles_for_tracking(particles)
        beamline.line.apply(lambda e: self.__generate_observation_points(e, beamline.length), axis=1)
        self.__add_input('run_track_beamline')
        self.__add_input('end_track')
        return self

    def __ptc_track(self, particles, beamline, **kwargs):
        if len(particles) == 0:
            print("No particles to track... Doing nothing.")
            return
        self.__add_input('ptc_create_universe')
        self.__add_input('ptc_create_layout', (False, 1, 4, 3, True))
        self.__add_particles_for_tracking(particles, True)
        beamline.line.apply(lambda e: self.__generate_observation_points_ptc(e, beamline.length), axis=1)
        self.__add_input('ptc_track', (
            5,
            0.0,
            False,
            True,
            1,
            True,
            True,
            'ptctrack',
            '.tfs'
        ))
        self.__add_input('ptc_track_end')
        self.__add_input('ptc_end')

    def __generate_observation_points(self, e, length):
        if not e['AT_EXIT'] == length:
            self.__add_input('observe', (e.name,))

    def __generate_observation_points_ptc(self, e, length):
        if not e['AT_EXIT'] == length and e['CLASS'] == 'MARKER':
            self.__add_input('ptc_observe', (e.name,))

    def show_beam(self):
        """Add a MAD-X `show beam` command."""
        self.__add_input('show_beam')
        return self

    def save_beta(self, **kwargs):
        """Add a MAD-X `save beta` command."""
        self.__add_input('save_beta', (kwargs.get("name", "BETA0"), kwargs.get("place", "#s")))
        return self

    def stop(self):
        """Add a MAD-X `stop` command (useful to act as a `break point`)."""
        self.__add_input('stop')
