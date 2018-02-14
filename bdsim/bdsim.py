import subprocess as sub
import re
import pandas as pd
import numpy as np
from ..simulator import Simulator
from ..lib import pybdsim


INPUT_FILENAME = 'input.bdsim'


class BdsimException(Exception):
    """Exception raised for errors in the Twiss module."""

    def __init__(self, m):
        self.message = m


class BDSim(Simulator):
    """A Python wrapper around the BDSim software.

    Uses `pybdsim` .
    """

    EXECUTABLE_NAME = 'bdsim'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._exec = BDSim.EXECUTABLE_NAME
        self._bdsim_machine = pybdsim.Builder.Machine()

    def run(self, **kwargs):
        """Run bdsim as a subprocess."""

        if self._get_exec() is None:
            raise BdsimException("Can't run BDSim if no valid path and executable are defined.")

        p = sub.Popen(f"{self._get_exec()} --file={INPUT_FILENAME}",
                      stdin=sub.PIPE,
                      stdout=sub.PIPE,
                      stderr=sub.STDOUT,
                      cwd=".",
                      shell=True
                      )
        self._output = p.communicate()[0].decode()
        self._warnings = [line for line in self._output.split('\n') if re.search('warning|fatal', line)]
        self._fatals = [line for line in self._output.split('\n') if re.search('fatal', line)]
        self._last_context = kwargs.get("context", {})
        if kwargs.get('debug', False):
            print(self._output)
        return self

    def beam_from_file(self, beam_file, **kwargs):
        b = pybdsim.Beam.Beam(
            particletype=kwargs.get('particletype', 'proton'),
            energy=kwargs['energy'],
            distrtype='userfile',
            distrFile=f"\"{beam_file}\"",
            distriFileFormat='"x[m]:xp[rad]:y[m]:yp[rad]:E[MeV]"'
        )
        self._bdsim_machine.AddBeam(b)