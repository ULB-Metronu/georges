import os
import re
import subprocess as sub

import jinja2
import numpy as np
import pandas as pd

from georges.simulator import Simulator
from .grammar import g4beamline_syntax

INPUT_FILENAME = 'input.g4bl'


class G4BeamlineException(Exception):
    """Exception raised for errors in the G4Beamline module."""

    def __init__(self, m):
        self.message = m


def element_to_g4beamline(e, fringe=None, build_solids=None, **kwargs):
    """Convert a pandas.Series representation onto a G4Beamline sequence element."""
    # For each element place a detector
    g4bl = "zntuple z={} file='Detector{}' format=ascii \n".format(e['AT_CENTER'] * 1000, e.name)
    if e['TYPE'] in ['QUADRUPOLE']:
        g4bl += g4beamline_syntax['quadrupole'].format(e.name,
                                                       e['LENGTH'] * 1000,
                                                       e['LENGTH'] * 1000,
                                                       e['APERTURE'] * 1000,
                                                       "{{{{ {} or '0.0' }}}}".format(e['CIRCUIT']))
        if not fringe:
            g4bl += "fringe=0"

        # if(kwargs.get("misalignment")):
        # else:
        g4bl += "\n place {} z={}\n".format(e.name, e['AT_CENTER'] * 1000)

    if e['TYPE'] in ['COLLIMATOR', 'SLITS']:
        if e['APERTYPE'] == 'CIRCLE':
            g4bl += g4beamline_syntax['ccoll'].format(e.name, e['APERTURE'] * 1000, e['LENGTH'] * 1000)
            g4bl += "\n place {} z={}\n".format(e.name, e['AT_CENTER'] * 1000)

        if e['APERTYPE'] == 'RECTANGLE':
            g4bl += g4beamline_syntax['rcoll'].format(e.name, e['LENGTH'] * 1000)

            g4bl += "\n place {} z={} {}={} rename={}_1\n".format(e.name, e['AT_CENTER'] * 1000,
                                                                  e['SLITS_PLANE'].lower(),
                                                                  f"(0.5*{{{{{e.name}_APERTURE or '0.1' }}}}*1000)+60/2",
                                                                  # do not forget the width of the slits
                                                                  e.name)

            g4bl += "\n place {} z={} {}={} rename={}_2\n".format(e.name, e['AT_CENTER'] * 1000,
                                                                  e['SLITS_PLANE'].lower(),
                                                                  f"(-0.5*{{{{{e.name}_APERTURE or '0.1' }}}}*1000)-60/2",
                                                                  # do not forget the width of the slits
                                                                  e.name)

    if e['TYPE'] == 'SOLIDS' and build_solids:

        if e['SOLIDS_FILE'] is None and e['SOLIDS_FILE'] is np.nan:
            raise G4BeamlineException(f"No files are provided for solids {e.name}.")

        solids_data = pd.read_csv(e['SOLIDS_FILE'])
        solidsinput = solids_data.apply(lambda b: g4beamline_syntax['tesselatedsolids'].format(b['NAME'],
                                                                                               b['MATERIAL'],
                                                                                               b['FILE']), axis=1)
        g4bl += '\n'.join(str(x) for x in solidsinput)

        soldisplacement = solids_data.apply(lambda b: g4beamline_syntax['place_solids'].format(b['NAME'],
                                                                                               b['POSX'],
                                                                                               b['POSY'] * 1000,
                                                                                               (b['POSZ'] + e[
                                                                                                   'AT_CENTER']) * 1000,
                                                                                               b['ROTX'],
                                                                                               b['ROTY'],
                                                                                               f"{{{{{e.name}_ROTATION or '0.0' }}}}"),
                                            axis=1)
        g4bl += '\n'
        g4bl += '\n'.join(str(x) for x in soldisplacement)

    return g4bl


def sequence_to_g4beamline(sequence, **kwargs):
    """Convert a pandas.DataFrame sequence onto a G4Beamline input."""
    sequence.sort_values(by='AT_CENTER', inplace=True)
    input = ""
    if sequence is None:
        return ""

    input += '\n'.join(sequence.apply(lambda e: element_to_g4beamline(e, **kwargs), axis=1)) + '\n'

    if sequence.iloc[-1]['TYPE'] == 'MARKER':  # If it is a marker : place a virtual detector to have results

        input += g4beamline_syntax['virtual_det'].format(sequence.iloc[-1]['AT_CENTER'] * 1000,
                                                         sequence.iloc[-1].name)
    return input


class G4Beamline(Simulator):
    """A Python wrapper around the G4Beamline executable.

    Sequence and command will be converted with the G4Beamline grammar and pipe'd to the subprocess.
    """

    EXECUTABLE_NAME = 'g4bl'

    def __init__(self, **kwargs):

        self._fringe = kwargs.get('fringe', False)
        self._build_degrader = kwargs.get('build_degrader', False)
        self.build_solids = kwargs.get('build_solids', None)

        super().__init__(**kwargs)

        self._exec = G4Beamline.EXECUTABLE_NAME

    def _attach(self, beamline):
        if beamline.length is None or pd.isnull(beamline.length):
            raise G4BeamlineException("Beamline length not defined.")

        self.__add_input('define_physics')
        self.__add_input('define_world')
        self.__add_input('keep_protons')
        self.__add_input('start_command')
        self._input += sequence_to_g4beamline(beamline.line, fringe=self._fringe,
                                              build_solids=self.build_solids)

    def _add__detector(self, e):
        self.__add_input('add_detector', (e['AT_CENTER'], e.name))

    def track(self, particles):
        if len(particles) == 0:
            print("No particles to track... Doing nothing.")
            return

        # Create the input file
        self.__add_input('beam_input', (len(particles),))
        with open('input_beam.dat', 'w+') as f:
            f.write('#BLTrackFile\n')
            f.write('# x y z Px Py Pz t PDGid EventID TrackID ParentID Weight\n')
            f.write('# mm mm mm MeV/c MeV/c MeV/c ns - - - - -\n')
            particles.to_csv(f, header=False, sep=' ', index=None
                             , columns=['X', 'Y', 'Z', 'PX', 'PY', 'PZ',
                                        't', 'PDGid', 'EventId', 'TrackId',
                                        'ParentId', 'Weight'])

    def __add_input(self, keyword, strings=()):
        self._input += g4beamline_syntax[keyword].format(*strings) + '\n'

    def run(self, **kwargs):
        """Run madx as a subprocess."""

        template_input = jinja2.Template(self._input).render(kwargs.get("context", {}))
        if kwargs.get("debug", False) >= 2:
            print(template_input)
        if self._exec is None:
            raise G4BeamlineException("Can't run G4Beamline if no valid path and executable are defined.")

        # Write the file for g4beamline
        file = open(INPUT_FILENAME, 'w')
        file.write(template_input)
        file.flush()
        file.close()

        cmd = " ".join([self._exec, INPUT_FILENAME])

        p = sub.Popen(cmd,
                      stdin=sub.PIPE,
                      stdout=sub.PIPE,
                      stderr=sub.STDOUT,
                      cwd=".",
                      shell=True
                      )
        self._output = p.communicate(input=template_input.encode())[0].decode()
        self._warnings = [line for line in self._output.split('\n') if re.search('Warning|Fatal Exception', line)]
        self._fatals = [line for line in self._output.split('\n') if re.search('Fatal Exception', line)]
        self._last_context = kwargs.get("context", {})
        if kwargs.get('debug', False):
            print(self._output)
        os.remove(INPUT_FILENAME)
        return self
