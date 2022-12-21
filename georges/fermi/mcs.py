"""
TODO
"""
import numpy as _np

from .. import ureg as _ureg


class ScatteringModelType(type):
    @staticmethod
    def t(pv: float, p1v1: float, **kwargs) -> float:
        return 0.0


class FermiRossi(metaclass=ScatteringModelType):
    """ """

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
        chi_0 = kwargs["material"].radiation_length.m_as("cm")
        return (es / pv) ** 2 * (1 / chi_0)


class DifferentialHighland(metaclass=ScatteringModelType):
    """ """

    @staticmethod
    def length(x, chi0):
        """

        Args:
            x:
            chi0:

        Returns:

        """
        return x / chi0

    @staticmethod
    def f_dh(length):
        """

        Args:
            length:

        Returns:

        """
        return 0.970 * (1 + (_np.log(length) / 20.7)) * (1 + (_np.log(length) / 22.7))

    @staticmethod
    def t(pv: float, p1v1: float, **kwargs) -> float:
        """

        Args:
            pv:
            p1v1:
            **kwargs:

        Returns:

        """
        material = kwargs.get("material")
        es = 14.1  # MeV
        chi0 = material.radiation_length.m_as("cm")
        x = material.required_thickness(kinetic_energy_out=pv * _ureg.MeV, kinetic_energy_in=p1v1 * _ureg.MeV).m_as(
            "cm",
        )
        return DifferentialHighland.f_dh(DifferentialHighland.length(x, chi0)) * (es / pv) ** 2 * (1 / chi0)


class ICRU(metaclass=ScatteringModelType):
    """ """

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
    """ """

    @staticmethod
    def t(pv: float, p1v1: float, **kwargs) -> float:
        """

        Args:
            pv:
            p1v1:
            **kwargs:

        Returns:

        """
        material = kwargs["material"]
        es = 15.0  # MeV
        chi_s = material.scattering_length.m_as("cm")
        return (es / pv) ** 2 * (1 / chi_s)


class DifferentialMoliere(metaclass=ScatteringModelType):
    """ """

    @staticmethod
    def t(pv: float, p1v1: float, **kwargs) -> float:
        """

        Args:
            pv:
            p1v1:
            **kwargs:

        Returns:

        """
        material = kwargs["material"]
        es = 15.0  # MeV
        chi_s = material.scattering_length.m_as("cm")
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
        return (
            0.5244
            + 0.1975 * _np.log10(1 - (pv / p1v1) ** 2)
            + 0.2320 * _np.log10(pv)
            - 0.0098 * _np.log10(pv) * _np.log10(1 - (pv / p1v1) ** 2)
        )
