import shutil
import subprocess as sub
import jinja2
import re
import pandas as pd
from madx.madx.grammar import madx_syntax

SUPPORTED_PROPERTIES = ['ANGLE', 'APERTYPE', 'E1', 'E2', 'FINT', 'HGAP', 'THICK', 'TILT']


def element_to_mad(e):
    """Convert a pandas.Series representation onto a MAD-X sequence element."""
    if not e['PHYSICAL'] or pd.isnull(e['PHYSICAL']):
        return ""
    mad = "{}: {}, ".format(e.name, e.CLASS)
    mad += ', '.join(["{}={}".format(p, e[p]) for p in SUPPORTED_PROPERTIES if pd.notnull(e[p])])
    if pd.notnull(e['ORBIT_LENGTH']): mad += ", L={}".format(e['ORBIT_LENGTH'])
    if pd.notnull(e['APERTYPE']): mad += ", APERTURE={}".format(str(e['APERTURE']).strip('[]'))
    if pd.notnull(e.get('PLUG')) and pd.notnull(e.get('CIRCUIT')): mad += ", {}:={}".format(e['PLUG'], e['CIRCUIT'])
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
        input += '\n'.join(sequence['CIRCUIT'].dropna().map(lambda c: "{}:={{{{ {} }}}};".format(c, c)))
    return input


def read_madx_twiss(file):
    """Read a MAD-X Twiss TFS file to a dataframe."""
    headers = pd.read_csv(file, skiprows=45, nrows=0, delim_whitespace=True)
    headers.drop(headers.columns[[0,1]], inplace=True, axis=1)
    df = pd.read_csv(file, header=None, names=headers, na_filter=False, skiprows=47, delim_whitespace=True)
    df.index.name = 'NAME'
    return df


def read_ptc_twiss(file):
    """Read a MAD-X PTC Twiss TFS file to a dataframe."""
    headers = pd.read_csv(file, skiprows=88, nrows=0, delim_whitespace=True)
    headers.drop(headers.columns[[0,1]], inplace=True, axis=1)
    df = pd.read_csv(file, header=None, names=headers, na_filter=False, skiprows=90, delim_whitespace=True)
    df.index.name = 'NAME'
    return df


def read_madx_tracking(file):
    """Read a MAD-X Tracking onetable=true file to a dataframe."""
    column_names = ['ID','TURN','X','PX','Y','PY','T','PT','S','E']
    data = pd.read_csv(file, skiprows=54, delim_whitespace=True, names=column_names)
    return data.apply(pd.to_numeric, errors="ignore").dropna()


def read_ptc_tracking(file):
    """Read a PTC Tracking 'one' file to a dataframe."""
    column_names = ['ID', 'TURN', 'X', 'PX', 'Y', 'PY', 'T', 'PT', 'S', 'E']
    data = pd.read_csv(file, skiprows=9, delim_whitespace=True,
                       names=column_names) \
              .apply(pd.to_numeric, errors="ignore").dropna()
    return data[data['TURN'] == 1]

class MadxException(Exception):
    """Exception raised for errors in the Madx module."""

    def __init__(self, m):
        self.message = m


class Madx:
    """A Python wrapper around the MAD-X executable.

    Sequence and command will be converted with the MAD-X grammar and pipe'd to the subprocess.
    """
    def __init__(self, **kwargs):
        self.__beamline = kwargs.get('beamline', None)
        if self.__beamline:
            self.__input = sequence_to_mad(self.__beamline.line)
        self._path = kwargs.get('path', "")
        self.__madx = kwargs.get('madx', None)
        self.__warnings = []
        self.__output = ""
        self.__template_input = ""

    def __get_madx_path(self):
        return self.__madx if self.__madx is not None else shutil.which("madx")

    def __add_input(self, keyword, strings=()):
        self.__input += madx_syntax[keyword].format(*strings) + '\n'

    def run(self, context):
        """Run madx as a subprocess."""
        self.__input += madx_syntax['stop']
        self.__template_input = jinja2.Template(self.__input).render(context)
        p = sub.Popen([self.__get_madx_path()],
                      stdin=sub.PIPE,
                      stdout=sub.PIPE,
                      stderr=sub.STDOUT,
                      cwd=self._path,
                      shell=True
                      )
        self.__output = p.communicate(input=self.__template_input.encode())[0].decode()
        self.__warnings = [line for line in self.__output.split('\n') if re.search('warning|fatal', line)]
        self.__fatals = [line for line in self.__output.split('\n') if re.search('fatal', line)]
        return self

    def print_input(self, context):
        """Print the rendered MAD-X input."""
        print(jinja2.Template(self.__input).render(context))

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
        return self.__template_input

    @property
    def output(self):
        """Return the output of the last MAD-X run."""
        return self.__output

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

    def beam(self):
        """Add a MAD-X `beam` command."""
        self.__add_input('beam')
        self.use_sequence(self.__beamline.name)
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

    def makethin(self, sequence, style='TEAPOT', dipole_slices=1, quadrupole_slices=3):
        """Add a MAD-X `makethin` command."""
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

    def track(self, particles, **kwargs):
        """Add a (ptc) `track` command."""
        if kwargs.get('ptc'):
            self.__ptc_track(particles, **kwargs)
        else:
            self.__track(particles)

    def __track(self, particles):
        if self.__beamline is None:
            print("No lattice defined.")
            return
        if len(particles) == 0:
            print("No particles to track... Doing nothing.")
            return
        self.makethin(self.__beamline.name)
        self.__add_input('track_beamline')
        self.__add_particles_for_tracking(particles)
        self.__beamline.line.apply(self.__generate_observation_points, axis=1)
        self.__add_input('run_track_beamline')
        self.__add_input('end_track')
        return self

    def __ptc_track(self, particles, **kwargs):
        if self.__beamline is None:
            print("No lattice defined.")
            return
        if len(particles) == 0:
            print("No particles to track... Doing nothing.")
            return
        self.__add_input('ptc_create_universe')
        self.__add_input('ptc_create_layout', (False, 1, 4, 3, True))
        self.__add_particles_for_tracking(particles, True)
        self.__beamline.line.apply(self.__generate_observation_points_ptc, axis=1)
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

    def __generate_observation_points(self, e):
        if not e['AT_EXIT'] == self.__beamline.length:
            self.__add_input('observe', (e.name,))

    def __generate_observation_points_ptc(self, e):
        if not e['AT_EXIT'] == self.__beamline.length and e['CLASS'] == 'MARKER':
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
