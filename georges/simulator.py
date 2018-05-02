import shutil
import jinja2
import os
import sys


class SimulatorException(Exception):
    """Exception raised for errors in the Simulator module."""

    def __init__(self, m):
        self.message = m


class Simulator:
    """Base class to support external 'simulator' programs (MAD-X, BDSim, etc.)."""

    def __init__(self, beamlines=None, path=None, **kwargs):
        """

        :param beamlines:
        :param path: path to the simulator executable (default: lookup using 'which')
        :param kwargs:
        """
        self._input = ""
        self._exec = ""
        self._output = None
        self._warnings = []
        self._fatals = []
        self._context = {}
        self._last_context = None
        self._path = path
        self._beamlines = []
        if beamlines and not isinstance(beamlines, list):
            raise SimulatorException("The 'beamlines' argument must be a list (if defined).")
        for b in beamlines:
            self._attach(b)

    def _attach(self, beamline):
        """Attach a beamline to the simulator instance."""
        print(beamline)
        self._beamlines.append(beamline)

    def _get_exec(self):
        """Retrive the path to the simulator executable."""
        if self._path is not None:
            return os.path.join(self._path, self._exec)
        elif: 'win' in sys.platform:
            return shutil.which(self._exec)
        else:
            return shutil.which(self._exec, path=".:/usr/local/bin:/usr/bin:/bin")

    def _add_input(self, keyword, *args, **kwargs):
        """Uses the simulator's own syntax to add to the input"""
        self._input += self._syntax[keyword].format(*args, **kwargs) + '\n'

    def print_warnings(self):
        """Print warnings from the previous execution run."""
        [print(w) for w in self._warnings]

    @property
    def warnings(self):
        """Return warnings from the previous execution run."""
        return self._warnings

    @property
    def fatals(self):
        """Return fatal errors from the previous execution run."""
        return self._fatals

    @property
    def input(self):
        """Return the current state of the input to be generated."""
        if self._last_context is None:
            raise SimulatorException("Asking for input but no context yet.")
        return jinja2.Template(self._input).render(self._last_context)

    @property
    def raw_input(self):
        """Return the current raw input."""
        return self._input

    @property
    def output(self):
        """Return the output of the last simulator run."""
        return self._output

    @property
    def beamlines(self):
        """Return the list of beamlines attached to the simulator."""
        return self._beamlines
