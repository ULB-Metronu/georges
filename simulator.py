import jinja2


class Simulator:
    @property
    def path(self):
        """Current simulator path."""
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
    def input(self, **kwargs):
        """Return the current MAD-X input."""
        rendered = kwargs.get("rendered", False)
        if rendered:
            return jinja2.Template(self.__input).render(self.context)
        else:
            return self.__input

    @property
    def output(self):
        """Return the output of the last simulator run."""
        return self.__output

    @property
    def beamlines(self):
        """Return the list of beamlines attached to the simulator."""
        return self.__beamlines

    @property
    def context(self):
        """The current state of the beamline."""
        return self.__context

    @context.setter
    def context(self, c):
        self.__context = c
