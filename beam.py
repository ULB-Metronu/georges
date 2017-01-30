import pandas as pd
import pandas.core.common
import numpy as np
import georges.physics as physics

PARTICLE_TYPES = {'proton', 'antiproton', 'electron', 'position'}
PHASE_SPACE_DIMENSIONS = ['X', 'PX', 'Y', 'PY', 'DPP', 'DT']


class BeamException(Exception):
    """Exception raised for errors in the Beam module."""

    def __init__(self, m):
        self.message = m


class Beam:
    """Particle beam to be tracked in a beamline or accelerator model.

    The internal representation is essentially a pandas DataFrame.
    """

    def __init__(self, *args, **kwargs):
        """
        :param args: distribution of particles to initialize the beam with. Should be pandas.DataFrame() friendly.
        :param kwargs: optional parameters include:
            - particle: identifier for the particle type (default: proton).
            - energy:
        """
        self.__distribution = None
        if len(args) >= 1:
            self.__initialize_distribution(args[0])

        self.__particle = None
        self.particle = kwargs.get('particle', 'proton')
        self.__energy = kwargs.get('energy', None)

    @property
    def distribution(self):
        """Return a dataframe containing the beam's particles distribution."""
        return self.__distribution

    @property
    def particle(self):
        """Return the particle type."""
        return self.__particle

    @particle.setter
    def particle(self, p):
        if p in PARTICLE_TYPES:
            self.__particle = p
        else:
            raise BeamException("Invalid particle type")

    @property
    def dims(self):
        """Return the dimensions of the beam's phase-space."""
        return self.__dims

    @property
    def n_particles(self):
        """Return the number of particles in the beam's distribution."""
        return self.__n_particles

    @property
    def energy(self):
        """Return the beam's energy in MeV."""
        if self.__energy is None:
            raise BeamException("Energy of the beam has not been set!")
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
        return self.__distribution.dropna().std()

    @property
    def halo(self):
        """Return a dataframe containing the 1st, 5th, 95th and 99th percentiles of each dimensions."""
        return pd.concat([
            self.__distribution.quantile(0.01),
            self.__distribution.quantile(0.05),
            self.__distribution.quantile(0.95),
            self.__distribution.quantile(0.99)
        ], axis=1).rename(columns={0: '1%', 1: '5%', 2: '95%', 3: '99%'})

    @property
    def coupling(self):
        """Return a dataframe containing the covariances (coupling) between each dimensions."""
        return self.__distribution.cov()

    def __getitem__(self, item):
        if item not in PHASE_SPACE_DIMENSIONS[:self.__dims]:
            raise BeamException("Trying to access an invalid data from a beam.")
        return self.__distribution[item]

    def __set_columns_names(self):
        self.__distribution.columns = PHASE_SPACE_DIMENSIONS[:self.__dims]

    def __initialize_distribution(self, distribution):
        """Try setting the internal pandas.DataFrame with a distribution."""
        try:
            self.__distribution = pd.DataFrame(distribution)
        except pd.core.common.PandasError:
            raise BeamException("Trying to initialize a beam from an invalid type.")
        self.__n_particles = self.__distribution.shape[0]
        self.__dims = self.__distribution.shape[1]
        if self.__dims < 2 or self.__dims > 6:
            raise BeamException("Trying to initialize a beam distribution with invalid dimensions.")
        self.__set_columns_names()

    def from_5d_multigaussian_distribution(self, n, **kwargs):
        """Initialize a beam with a 5D particle distribution."""
        keys = {'X', 'PX', 'Y', 'PY', 'DPP', 'XRMS', 'PXRMS', 'YRMS', 'PYRMS', 'DPPRMS'}
        if any([k not in keys for k in kwargs.keys()]):
            raise BeamException("Invalid argument for a multigaussian distribution.")
        self.__initialize_distribution(pd.DataFrame(np.random.multivariate_normal(
            [kwargs.get('X', 0),
             kwargs.get('PX', 0),
             kwargs.get('Y', 0),
             kwargs.get('PY', 0),
             kwargs.get('DPP', 0)
             ],
            np.diag([kwargs.get('XRMS', 0),
                     kwargs.get('PXRMS', 0),
                     kwargs.get('YRMS', 0),
                     kwargs.get('PYRMS', 0),
                     kwargs.get('DPPRMS', 0)
                     ]) ** 2,
            n
        )))
        self.__set_columns_names()
        return self
