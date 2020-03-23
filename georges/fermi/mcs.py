"""
TODO
"""
import numpy as _np


def radiation_length(material) -> float:
    """

    Args:
        material:

    Returns:

    """
    a: float = material.atomic_a
    z: float = material.atomic_z
    return 716.4 * a * (1 / (z * (z + 1) * (_np.log(287 / _np.sqrt(z)))))


def scattering_length(material) -> float:
    """
    See "Techniques of Proton Radiotherapy:Transport Theory", B. Gottschalk, 2012.

    Args:
        material:

    Returns:

    """
    rho = material.density
    if material is 'water':
        return 46.88/rho
    if material is 'air':
        return 46.76 / rho
    alpha: float = 0.0072973525664  # Fine structure constant
    avogadro: float = 6.02e23  # Avogadro's number
    re: float = 2.817940e-15 * 100  # Classical electron radius (in cm)
    a: float = material.atomic_a
    z: float = material.atomic_z
    return 1 / (rho * alpha * avogadro * re ** 2 * z ** 2 * (2 * _np.log(33219 * (a * z) ** (-1 / 3)) - 1) / a)


class ScatteringModelType(type):
    @staticmethod
    def t(pv: float, p1v1: float, **kwargs) -> float:
        return 0.0


class FermiRossi(metaclass=ScatteringModelType):
    """

    """
    @staticmethod
    def t(pv: float, p1v1: float, **kwargs) -> float:
        """

        Args:
            pv:
            p1v1:
            **kwargs:

        Returns:

        """
        es = 15.0  # MeV
        chi_0 = 19.32
        return (es/pv) ** 2 * (1/chi_0)


class DifferentialHighland(metaclass=ScatteringModelType):
    """

    """

    @staticmethod
    def l(x, chi0):
        """

        Args:
            x:
            chi0:

        Returns:

        """
        return x / chi0

    @staticmethod
    def f_dh(l):
        """

        Args:
            l:

        Returns:

        """
        return 0.970 * (1 + (_np.log(l) / 20.7)) * (1 + (_np.log(l) / 22.7))

    @staticmethod
    def t(pv: float, p1v1: float, **kwargs) -> float:
        """

        Args:
            pv:
            p1v1:
            **kwargs:

        Returns:

        """
        material = kwargs.get('material')
        es = 14.1  # MeV
        chi0 = radiation_length(material=material)
        x = material.thickness(k_out=pv, k_in=p1v1)
        return DifferentialHighland.f_dh(DifferentialHighland.l(x, chi0)) * (es / pv) ** 2 * (1 / chi0)


class ICRU(metaclass=ScatteringModelType):
    """

    """
    @staticmethod
    def t(pv: float, p1v1: float, **kwargs) -> float:
        """

        Args:
            pv:
            p1v1:
            **kwargs:

        Returns:

        """
        pass


class ICRUProtons(metaclass=ScatteringModelType):
    """

    """
    @staticmethod
    def t(pv: float, p1v1: float, **kwargs) -> float:
        """

        Args:
            pv:
            p1v1:
            **kwargs:

        Returns:

        """
        material = kwargs['material']
        es = 15.0  # MeV
        chi_s = scattering_length(material=material)
        return (es / pv) ** 2 * (1 / chi_s)


class DifferentialMoliere(metaclass=ScatteringModelType):
    """

    """
    @staticmethod
    def t(pv: float, p1v1: float, **kwargs) -> float:
        """

        Args:
            pv:
            p1v1:
            **kwargs:

        Returns:

        """
        material = kwargs['material']
        es = 15.0  # MeV
        chi_s = scattering_length(material=material)
        return DifferentialMoliere.f_dm(p1v1, pv) * (es / pv) ** 2 * (1 / chi_s)

    @staticmethod
    def f_dm(p1v1: float, pv: float):
        """

        Args:
            p1v1:
            pv:

        Returns:

        """
        if pv <= 0:
            raise ValueError("'pv' must be > 0.")
        if p1v1 <= 0:
            raise ValueError("'p1v1' must be > 0.")
        if p1v1 <= pv:
            raise ValueError("Initial 'p1v1' must be larger than final 'pv'.")
        return 0.5244 \
               + 0.1975 * _np.log10(1 - (pv / p1v1) ** 2) \
               + 0.2320 * _np.log10(pv) \
               - 0.0098 * _np.log10(pv) * _np.log10(1 - (pv / p1v1) ** 2)
