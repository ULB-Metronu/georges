import numpy as np


def inverse_scattering_length(**kwargs):
    db = kwargs.get('db')
    material = kwargs['material']
    return 1.0/0.03372681281618887


class FermiRossi:
    """"""
    @staticmethod
    def t(pv, p1v1, **kwargs):
        es = 15.0  # MeV
        chi_0 = 19.32
        return (es/pv) ** 2 * (1/chi_0)


class DifferentialMoliere:
    """"""
    @staticmethod
    def t(pv, p1v1, **kwargs):
        db = kwargs.get('db')
        material = kwargs['material']
        es = 15.0  # MeV
        chi_s = inverse_scattering_length(material=material, db=db)
        return DifferentialMoliere.f_dm(p1v1, pv) * (es / pv) ** 2 * (1 / chi_s)

    @staticmethod
    def f_dm(p1v1, pv):
        return 0.5244 \
               + 0.1975 * np.log10(1 - (pv / p1v1) ** 2) \
               + 0.2320 * np.log10(pv) \
               - 0.0098 * np.log10(pv) * np.log10(1 - (pv / p1v1) ** 2)
