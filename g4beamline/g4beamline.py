import shutil
import subprocess as sub
import jinja2
import re
import os
import pandas as pd
from ..simulator import Simulator
from .grammar import g4beamline_syntax


class G4BeamlineException(Exception):
    """Exception raised for errors in the G4Beamline module."""

    def __init__(self, m):
        self.message = m


def element_to_g4beamline(e):
    """Convert a pandas.Series representation onto a G4Beamline sequence element."""
    g4bl = ""

    # For each element place a detector
    g4bl = "place detector z={} rename=Detector{}\n".format(e['AT_CENTER'], e.name)
    if e['TYPE'] in ['QUADRUPOLE']:
        g4bl+="genericquad {} " \
              "fieldLength={} " \
              "ironLength={} " \
              "ironRadius=238 " \
              "apertureRadius={} " \
              "gradient={}*$Brho " \
              "ironMaterial=Fe " \
              "fieldMaterial=Vacuum " \
              "ironColor=1,0,0 " \
              "kill=1 " \
              "fringe={} " \
              "openAperture=1\n".format(e.name,
                                      e['LENGTH']*1000,
                                      e['LENGTH']*1000,
                                      e['APERTURE']*1000,
                                      "{{{{ {} or '0.0' }}}}".format(e['CIRCUIT']), ## To change for generic case
                                      0) ##kwargs.get("fringe")

        #if(kwargs.get("misalignment")):
        g4bl+="place {} z={}\n".format(e.name,e['AT_CENTER']*1000)

    return g4bl

def sequence_to_g4beamline(sequence):
    """Convert a pandas.DataFrame sequence onto a G4Beamline input."""
    sequence.sort_values(by='AT_CENTER', inplace=True)
    input=""
    if sequence is None:
        return ""

    input+= '\n'.join(sequence.apply(element_to_g4beamline, axis=1)) + '\n'

    return input


class G4Beamline(Simulator):
    """A Python wrapper around the G4Beamline executable.

    Sequence and command will be converted with the G4Beamline grammar and pipe'd to the subprocess.
    """

    EXECUTABLE_NAME = 'g4bl'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _attach(self, beamline):
        if beamline.length is None or pd.isnull(beamline.length):
            raise G4BeamlineException("Beamline length not defined.")

        self.__add_input('define_Brho')
        self.__add_input('define_detector')
        self._input += sequence_to_g4beamline(beamline.line)

    def _add__detector(self,e):
        self.__add_input('add_detector', (e['AT_CENTER'],e.name))

    def __add_particles_for_tracking(self,particles):
        if {'X', 'PX', 'Y', 'PY', 'DPP'} > set(particles):
            return
        for r in particles.iterrows():
            self.__add_input('beam_start', tuple(r[1]))

    def track(self, particles):
        if len(particles) == 0:
            print("No particles to track... Doing nothing.")
            return

        self.__add_particles_for_tracking(particles)

    def __add_input(self, keyword, strings=()):
        self._input += g4beamline_syntax[keyword].format(*strings) + '\n'

    def run(self, **kwargs):
        """Run madx as a subprocess."""

        template_input = jinja2.Template(self._input).render(kwargs.get("context", {}))
        if kwargs.get("debug", False) >= 2:
            print(template_input)
        if self._get_exec() is None:
            raise ("Can't run G4Beamline if no valid path and executable are defined.")

        # Write the file for g4beamline
        file = open('input_g4beamline.g4bl', 'w')
        file.write(template_input)
        file.flush()
        file.close()
        g4cmd=self._get_exec()+' input_g4beamline.g4bl'
        p = sub.Popen([g4cmd],
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
        #os.remove('input_g4beamline.g4bl')
        return self


