import subprocess as sub
import jinja2
import re
import pandas as pd
import numpy as np
from ..simulator import Simulator
from .grammar import bdsim_syntax

SUPPORTED_CLASSES = ['DRIFT',
                     'RBEND',
                     'SBEND',
                     'QUADRUPOLE',
                     'SEXTUPOLE',
                     'OCTUPOLE',
                     'DECAPOLE',
                     'MULTIPOLE',
                     'THINMULTIPOLE',
                     'VKICKER',
                     'HKICKER',
                     'KICKER',
                     'RF',
                     'RCOL',
                     'ECOL',
                     'DEGRADER',
                     'MUSPOILER',
                     'SHIELD',
                     'SOLENOID',
                     'LASER',
                     'TRANSFORM3D',
                     'ELEMENT',
                     'MARKER'
                     ]
SUPPORTED_PROPERTIES = ['ANGLE', 'APERTYPE', 'E1', 'E2', 'FINT', 'HGAP', 'THICK', 'TILT']
INPUT_FILENAME = 'input.bdsim'


class BdsimException(Exception):
    """Exception raised for errors in the Twiss module."""

    def __init__(self, m):
        self.message = m


def split_rbends(line, n=20):
    split_line = pd.DataFrame()
    for index, row in line.iterrows():
        if row['CLASS'] == 'RBEND' and pd.isnull(row.get('SPLIT')):
            angle = row['ANGLE'] / n
            length = row['L'] / n
            for i in range(0,n):
                row = row.copy()
                row.name = index + "_{}".format(i)
                row['SPLIT'] = True
                row['ANGLE'] = angle
                row['L'] = length
                split_line = split_line.append(row)
        else:
            split_line = split_line.append(row)
    split_line[['THICK']] = split_line[['THICK']].applymap(bool)
    return split_line


def element_to_bdsim(e):
    """Convert a pandas.Series representation onto a BDSim sequence element."""
    bdsim = ""
    if e.KEYWORD in ['MARKER', 'INSTRUMENT']:
        bdsim = "{}: {};".format(e.name.replace('$', ''), "marker")
    if e.KEYWORD in ['DRIFT', 'QUADRUPOLE', 'RBEND', 'SBEND']:
        bdsim = "{}: {}, l={}*m".format(e.name.replace('$', ''), e.KEYWORD.lower(), e.L)
        if e.get('BENDING_ANGLE') is not None and not np.isnan(e['BENDING_ANGLE']):
            bdsim += f",angle=-{e['BENDING_ANGLE']}"
        elif e.get('ANGLE') is not None and not np.isnan(e['ANGLE']):
            bdsim += f",angle=-{e.get('ANGLE', 0)}"
        else:
            # Angle property not supported by the element or absent
            bdsim += ""
        #if pd.notnull(e['APERTYPE']):
        #    bdsim += ", aperture={}*m".format(str(e['APERTURE']).strip('[]'))
        if pd.notnull(e.get('PLUG')) and pd.notnull(e.get('CIRCUIT')):
            bdsim += ", {}={{{{ {} or '0.0' }}}}".format(e['PLUG'].lower(), e['CIRCUIT'])
        bdsim += ';'
    return bdsim


def sequence_to_bdsim(seq, **kwargs):
    """Convert a pandas.DataFrame sequence onto a BDSim input."""
    if seq is None:
        return ""
    sequence = seq.line
    sequence.sort_values(by='S', inplace=True)
    if kwargs.get("split_rbend", False):
        sequence = split_rbends(sequence)
    # Drift smaller than 5 microns are discarded
    sequence = sequence.query("L > 5e-6 or KEYWORD != 'DRIFT'")
    i = "\n".join(
        sequence.reset_index()
                .drop_duplicates(subset='index', keep='last')
                .set_index('index')
                .apply(
            element_to_bdsim, axis=1)
    )

    i += "\n{}: line = ({});\n".format(seq.name, ",".join(sequence.index.map(lambda x: x.replace('$', ''))))
    return i


class BDSim(Simulator):
    """A Python wrapper around the BDSim software.

    Sequence and command will be converted with the BDSim grammar and pipe'd to the subprocess.
    """

    EXECUTABLE_NAME = 'bdsim'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._syntax = bdsim_syntax
        self._exec = BDSim.EXECUTABLE_NAME

    def attach(self, beamline, *args, **kwargs):
        super()._attach(beamline)
        self._input = "BRHO=2.3114; DEGREE=pi/180.0;"

        self._input += sequence_to_bdsim(beamline, split_rbend=False)
        self._add_input("use", line=beamline.name)

        for b in args:
            b.line.index = b.name + b.line.index
            self._input += sequence_to_bdsim(b, split_rbend=False)
            self._add_input("placement",
                            line=b.name,
                            reference_element=kwargs.get("placement").get(b.name),
                            reference_element_number=0,
                            x_placement=0,
                            z_placement=1
                            )

    def run(self, **kwargs):
        """Run bdsim as a subprocess."""
        self._add_input("options",
                        beampiperadius=32.5,
                        aperturetype="circular",
                        beampipethickness=1.0,
                        beampipematerial="Aluminium"
                        )
        self._add_input("beam",
                        particle='proton',
                        energy=230 + 938.272,
                        )

        template_input = jinja2.Template(self._input).render(kwargs.get('context', {}))
        if self._get_exec() is None:
            raise BdsimException("Can't run BDSim if no valid path and executable are defined.")
        with open(INPUT_FILENAME, 'w') as f:
            f.write(template_input)
        p = sub.Popen(f"{self._get_exec()} --file={INPUT_FILENAME}",
                      stdin=sub.PIPE,
                      stdout=sub.PIPE,
                      stderr=sub.STDOUT,
                      cwd=".",
                      shell=True
                      )
        self._output = p.communicate(input=template_input.encode())[0].decode()
        self._warnings = [line for line in self._output.split('\n') if re.search('warning|fatal', line)]
        self._fatals = [line for line in self._output.split('\n') if re.search('fatal', line)]
        self._last_context = kwargs.get("context", {})
        if kwargs.get('debug', False):
            print(self._output)
        #os.remove(INPUT_FILENAME)
        return self
