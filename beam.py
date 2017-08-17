import pandas as pd
import pandas.core.common
import numpy as np
from . import physics

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

        self.__particle = kwargs.get('particle', 'proton')
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
        return self.__distribution.std()

    @property
    def emit(self):
        """Return the emittance of the beam in both planes"""
        return {'X': np.sqrt(np.linalg.det(self.__distribution.head(len(self.__distribution))[['X', 'PX']].cov())),
                'Y': np.sqrt(np.linalg.det(self.__distribution.head(len(self.__distribution))[['Y', 'PY']].cov()))
                }

    @property
    def sigma(self):
        """Return the sigma matrix of the beam"""
        return self.__distribution.cov()

    @property
    def halo(self):
        """Return a dataframe containing the 1st, 5th, 95th and 99th percentiles of each dimensions."""

        return pd.concat([
            self.__distribution.quantile(0.01),
            self.__distribution.quantile(0.05),
            self.__distribution.quantile(1.0-0.842701),
            self.__distribution.quantile(0.842701),
            self.__distribution.quantile(0.95),
            self.__distribution.quantile(0.99)
        ], axis=1).rename(columns={0.01: '1%', 0.05: '5%', 1.0-0.842701: '20%', 0.842701: '80%', 0.95: '95%', 0.99: '99%'})

    @property
    def coupling(self):
        """Return a dataframe containing the covariances (coupling) between each dimensions."""
        return self.__distribution.cov()

    def __getitem__(self, item):
        if item not in PHASE_SPACE_DIMENSIONS[:self.__dims]:
            raise BeamException("Trying to access an invalid data from a beam.")
        return self.__distribution[item]

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
        self.__distribution.columns = PHASE_SPACE_DIMENSIONS[:self.__dims]

    def from_5d_multigaussian_distribution(self, n, **kwargs):
        """Initialize a beam with a 5D particle distribution."""
        keys = {'X', 'PX', 'Y', 'PY', 'DPP', 'XRMS', 'PXRMS', 'YRMS', 'PYRMS', 'DPPRMS'}
        if any([k not in keys for k in kwargs.keys()]):
            raise BeamException("Invalid argument for a multigaussian distribution.")
        self.from_5d_sigma_matrix(n,
                                  X=kwargs.get('X', 0),
                                  PX=kwargs.get('PX', 0),
                                  Y=kwargs.get('Y', 0),
                                  PY=kwargs.get('PY', 0),
                                  DPP=kwargs.get('DPP', 0),
                                  DPPRMS=kwargs.get('DPPRMS', 0),
                                  s11=kwargs.get('XRMS', 0),
                                  s12=0,
                                  s22=kwargs.get('PXRMS', 0),
                                  s33=kwargs.get('YRMS', 0),
                                  s34=0,
                                  s44=kwargs.get('PYRMS', 0)
                                  )
        return self

    def from_5d_sigma_matrix(self, n, **kwargs):
        """Initialize a beam with a 5D particle distribution from a \Sigma matrix."""
        s11 = kwargs.get('s11', 0)
        s12 = kwargs.get('s12', 0)
        s21 = s12
        s22 = kwargs.get('s22', 0)
        s33 = kwargs.get('s33', 0)
        s34 = kwargs.get('s34', 0)
        s43 = s34
        s44 = kwargs.get('s44', 0)
        sdpp = kwargs.get('DPPRMS', 0)
        self.__initialize_distribution(pd.DataFrame(np.random.multivariate_normal(
            [kwargs.get('X', 0),
             kwargs.get('PX', 0),
             kwargs.get('Y', 0),
             kwargs.get('PY', 0),
             kwargs.get('DPP', 0)
             ],
            np.array([
                [s11, s12, 0, 0, 0],
                [s21, s22, 0, 0, 0],
                [0, 0, s33, s34, 0],
                [0, 0, s43, s44, 0],
                [0, 0, 0, 0, sdpp]
            ]),
            n
        )))
        self.__distribution.columns = PHASE_SPACE_DIMENSIONS[:self.__dims]
        return self
