import pandas as pd
import pandas.core.common
import numpy as np
#from . import physics
from scipy.optimize import curve_fit

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

        #self.__distribution = None
        self.__distribution = kwargs.get('distribution', None)
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
    def twiss(self):
        """Return the Twiss parameters of the beam"""
        s11 = self.sigma['X']['X']
        s12 = self.sigma['X']['PX']
        s22 = self.sigma['PX']['PX']
        s33 = self.sigma['Y']['Y']
        s34 = self.sigma['Y']['PY']
        s44 = self.sigma['PY']['PY']
        return {
            'beta_x': s11 / self.emit['X'],
            'alpha_x': -s12 / self.emit['X'],
            'gamma_x': s22 / self.emit['X'],
            'beta_y': s33 / self.emit['Y'],
            'alpha_y': -s34 / self.emit['Y'],
            'gamma_y': s44 / self.emit['Y'],
        }

    @property
    def std_bpm(self):
        return self._std_bpm()

    def _std_bpm(self, **kwargs):
        """TODO"""
        def gaussian(x, a, mu, sigma):
            return a * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)) / (np.sqrt(2*np.pi)*sigma)

        def fit_bpm(d, ax=None):
            bs = np.array(
                [-31, -19.8, -15.8, -11.8, -7.8, -5.8, -3.8, -1.8, 0.0, 1.8, 3.8, 5.8, 7.8, 11.8, 15.8, 19.8, 31]) / 1000
            bsp = (bs[1:] + bs[:-1]) / 2
            w = 1.0 / (bs[1:] - bs[:-1])
            w[0] *= 0.7
            w[-1] *= 0.7
            hist = np.histogram(d, bs)
            x = bsp
            y = w * hist[0]
            ar = np.trapz(y / np.sum(y) * len(y), x)
            mean = np.mean(x * y / np.sum(y) * len(y))
            rms = np.std(x * y / np.sum(y) * len(y))
            popt, pcov = curve_fit(gaussian, x, y,
                                   p0=[ar, mean, rms],
                                   maxfev=10000,
                                   bounds=(
                                       (-np.inf, mean - 0.1 * np.abs(mean), 0.5 * rms),
                                       (np.inf, mean + 0.1 * np.abs(mean), 2 * rms)
                                   )
                                   )

            if ax is not None:
                kwargs['ax'].plot(bsp, w * hist[0], '*-')
                kwargs['ax'].plot(x, gaussian(x, *popt), 'ro:', label='fit')

            return [
                np.abs(popt[2]),
                np.sqrt(pcov[2, 2]) if pcov[2, 2] > 0 else 0.0
            ]

        # Simulate orbit trajectory shifts of 1mm to get statistical error on beam size
        data_x = np.array([fit_bpm(self.__distribution['X'], kwargs.get('ax')),
                           fit_bpm(self.__distribution['X'] + 0.002),
                           fit_bpm(self.__distribution['X'] - 0.002)])
        mean_x = np.mean(data_x[:, 0])
        error_x = np.sqrt(np.std(data_x[:, 0])**2 + np.max(data_x[:, 1])**2)
        data_y = np.array([fit_bpm(self.__distribution['Y'], kwargs.get('ax')),
                           fit_bpm(self.__distribution['Y'] + 0.002),
                           fit_bpm(self.__distribution['Y'] - 0.002)])
        mean_y = np.mean(data_y[:, 0])
        error_y = np.sqrt(np.std(data_y[:, 0])**2 + np.max(data_y[:, 1])**2)

        return {'X': [mean_x, error_x], 'Y': [mean_y, error_y]}

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
        ], axis=1).rename(columns={0.01: '1%',
                                   0.05: '5%',
                                   1.0-0.842701: '20%',
                                   0.842701: '80%',
                                   0.95: '95%',
                                   0.99: '99%'
                                   }
                          )

    @property
    def coupling(self):
        """Return a dataframe containing the covariances (coupling) between each dimensions."""
        return self.__distribution.cov()

    def __getitem__(self, item):
        if item not in PHASE_SPACE_DIMENSIONS[:self.__dims]:
            raise BeamException("Trying to access an invalid data from a beam.")
        return self.__distribution[item]

    def from_file(self, filename):
        ##
        print('Open file : ', filename)
         #beam=

        return self

    def __initialize_distribution(self, distribution, **kwargs):
        """Try setting the internal pandas.DataFrame with a distribution."""
        try:
            self.__distribution = pd.DataFrame(distribution)
        except pd.core.common.PandasError:
            if kwargs.get("filename"):
                self._from_file(kwargs.get("filename"))
            else:
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
                                  s11=kwargs.get('XRMS', 0)**2,
                                  s12=0,
                                  s22=kwargs.get('PXRMS', 0)**2,
                                  s33=kwargs.get('YRMS', 0)**2,
                                  s34=0,
                                  s44=kwargs.get('PYRMS', 0)**2
                                  )
        return self

    def from_twiss_parameters(self, n, **kwargs):
        """Initialize a beam with a 5D particle distribution from Twiss parameters."""
        keys = {'X', 'PX', 'Y', 'PY', 'DPP', 'DPPRMS', 'BETAX', 'ALPHAX', 'BETAY', 'ALPHAY', 'EMITX', 'EMITY'}
        if any([k not in keys for k in kwargs.keys()]):
            raise BeamException("Invalid argument for a twiss distribution.")
        betax = kwargs.get('BETAX', 1)
        alphax = kwargs.get('ALPHAX', 0)
        gammax = (1+alphax**2)/betax
        betay = kwargs.get('BETAY', 1)
        alphay = kwargs.get('ALPHAY', 0)
        gammay = (1 + alphay ** 2) / betay

        self.from_5d_sigma_matrix(n,
                                  X=kwargs.get('X', 0),
                                  PX=kwargs.get('PX', 0),
                                  Y=kwargs.get('Y', 0),
                                  PY=kwargs.get('PY', 0),
                                  DPP=kwargs.get('DPP', 0),
                                  DPPRMS=kwargs.get('DPPRMS', 0),
                                  s11=betax * kwargs['EMITX'],
                                  s12=-alphax * kwargs['EMITX'],
                                  s22=gammax * kwargs['EMITX'],
                                  s33=betay * kwargs['EMITY'],
                                  s34=-alphay * kwargs['EMITY'],
                                  s44=gammay * kwargs['EMITY']
                                  )
        return self

    def from_5d_sigma_matrix(self, n, **kwargs):
        """Initialize a beam with a 5D particle distribution from a \Sigma matrix."""
        s11 = kwargs.get('s11', 0)
        s12 = kwargs.get('s12', 0)
        s13 = kwargs.get('s13', 0)
        s14 = kwargs.get('s14', 0)
        s15 = kwargs.get('s15', 0)
        s21 = s12
        s22 = kwargs.get('s22', 0)
        s23 = kwargs.get('s23', 0)
        s24 = kwargs.get('s24', 0)
        s25 = kwargs.get('s25', 0)
        s31 = s13
        s32 = s23
        s33 = kwargs.get('s33', 0)
        s34 = kwargs.get('s34', 0)
        s35 = kwargs.get('s35', 0)
        s41 = s14
        s42 = s24
        s43 = s34
        s44 = kwargs.get('s44', 0)
        s45 = kwargs.get('s45', 0)
        s51 = s15
        s52 = s25
        s53 = s35
        s54 = s45
        s55 = kwargs.get('DPPRMS', 0)**2

        self.__initialize_distribution(pd.DataFrame(np.random.multivariate_normal(
            [kwargs.get('X', 0),
             kwargs.get('PX', 0),
             kwargs.get('Y', 0),
             kwargs.get('PY', 0),
             kwargs.get('DPP', 0)
             ],
            np.array([
                [s11, s12, s13, s14, s15],
                [s21, s22, s23, s24, s25],
                [s31, s32, s33, s34, s35],
                [s41, s42, s43, s44, s45],
                [s51, s52, s53, s54, s55]
            ]),
            n
        )))
        self.__distribution.columns = PHASE_SPACE_DIMENSIONS[:self.__dims]
        return self
