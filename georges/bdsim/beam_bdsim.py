import pandas as pd
import numpy as np
from ..beam import Beam
from .. import physics


PHASE_SPACE_DIMENSIONS = ['X', 'PX', 'Y', 'PY', 'P', 'ENERGY', 'PARENT_ID', 'PDG_ID', 'WEIGHT']


class BeamBdsimException(Exception):
    """Exception raised for errors in the Madx module."""

    def __init__(self, m):
        self.message = m


class BeamBdsim(Beam):

    """Particle beam resulting of a bdsim simulation.

    The internal representation is essentially a pandas DataFrame.
    """

    def __init__(self, distribution=None, *args, **kwargs):
        """

        :param distribution: distribution of particles to initialize the beam with. Should be pandas.DataFrame() friendly.
        :param particle: the particle type (default: 'proton', must be 'proton', 'antiproton', 'electron' or 'positron').
        :param energy: the reference energy of the beam
        :param args: optional parameters.
        :param kwargs: optional keyword parameters.
        """
        try:
            self.__initialize_distribution(distribution, *args, **kwargs)
        except BeamBdsimException:
            self.__dims = 9
            self.__distribution = pd.DataFrame(np.zeros((1, 9)))
            self.__distribution.columns = PHASE_SPACE_DIMENSIONS[:self.__dims]

    @property
    def distribution(self):
        """Return a dataframe containing the beam's particles distribution."""
        return self.__distribution

    @property
    def n_particles(self):
        """Return the number of particles in the beam's distribution."""
        return self.__n_particles

    @property
    def energy(self):
        """Return the beam's energy in MeV."""
        if self.__energy is None:
            raise BeamBdsimException("Energy of the beam has not been set!")
        return self.__energy

    @property
    def pc(self):
        """Return the beam's momentum in MeV."""
        return physics.energy_to_momentum(self.energy)

    @property
    def betarel(self):
        """Return the beam's relativistic beta value."""
        return physics.energy_to_beta(self.energy)

    @property
    def mean(self):
        """Return a dataframe containing the first order moments of each dimensions."""
        return self.__distribution.mean()

    @property
    def std(self):
        """Return a dataframe containing the second order moments of each dimensions."""
        return self.__distribution.std()

    @property
    def emit(self):
        """Return the emittance of the beam in both planes"""
        return {
            'X': np.sqrt(np.linalg.det(self.__distribution.head(len(self.__distribution))[['X', 'PX']].cov())),
            'Y': np.sqrt(np.linalg.det(self.__distribution.head(len(self.__distribution))[['Y', 'PY']].cov()))
        }

    def __initialize_distribution(self, distribution=None, *args, **kwargs):
        """Try setting the internal pandas.DataFrame with a distribution."""
        if distribution is not None:
            self.__distribution = distribution
        else:
            try:
                self.__distribution = pd.DataFrame(args[0])
            except (IndexError, ValueError):
                if kwargs.get("filename") is not None:
                    self.__distribution = Beam.from_file(kwargs.get('filename'), path=kwargs.get('path', ''))
                else:
                    return
        self.__n_particles = self.__distribution.shape[0]
        if self.__n_particles <= 0:
            raise BeamBdsimException("Trying to initialize a beam distribution with invalid number of particles.")
        self.__dims = self.__distribution.shape[1]
        if self.__dims < 2 or self.__dims > 9:
            raise BeamBdsimException("Trying to initialize a beam distribution with invalid dimensions.")
        self.__distribution.columns = PHASE_SPACE_DIMENSIONS[:self.__dims]

    @property
    def halo(self):
        """Return a dataframe containing the 1st, 5th, 95th and 99th percentiles of each dimensions."""

        return pd.concat([
            self.__distribution.quantile(0.01),
            self.__distribution.quantile(0.05),
            self.__distribution.quantile(1.0 - 0.842701),
            self.__distribution.quantile(0.842701),
            self.__distribution.quantile(0.95),
            self.__distribution.quantile(0.99)
        ], axis=1).rename(columns={0.01: '1%',
                                   0.05: '5%',
                                   1.0 - 0.842701: '20%',
                                   0.842701: '80%',
                                   0.95: '95%',
                                   0.99: '99%'
                                   }
                          )

