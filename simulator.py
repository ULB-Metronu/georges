
import shutil
import subprocess as sub
import jinja2
import re
import pandas as pd
f
class Madx(Simulator):
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





    def set(self, k, v):
        """Set a single variable in the context. Allows method chaining."""
        self.__context[k] = v
        return self

    def print_warnings(self):
        """Print warnings from the previous execution run."""
        [print(w) for w in self.__warnings]



class Simulator:
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