import numpy as np
from .stopping import thickness


def radiation_length(db, material):
    # See "reference"
    a = db.a(material)
    z = db.z(material)
    return 716.4 * a * (1 / (z * (z + 1) * (np.log(287 / np.sqrt(z)))))


def scattering_length(db, material):
    # See "Techniques of Proton Radiotherapy:Transport Theory", B. Gottschalk, 2012.
    rho = db.density(material)
    if material is 'water':
        return 46.88/rho
    if material is 'air':
        return 46.76 / rho
    alpha = 0.0072973525664  # Fine structure constant
    avogadro = 6.02e23  # Avogadro's number
    re = 2.817940e-15 * 100  # Classical electron radius (in cm)
    a = db.a(material)
    z = db.z(material)
    return 1 / (rho * alpha * avogadro * re ** 2 * z ** 2 * (2 * np.log(33219 * (a * z) ** (-1 / 3)) - 1) / a)


class FermiRossi:
    """"""
    @staticmethod
    def t(pv: float, p1v1: float, **kwargs):
        es = 15.0  # MeV
        chi_0 = 19.32
        return (es/pv) ** 2 * (1/chi_0)


class DifferentialHighland:
    """"""

    @staticmethod
    def l(x, chi0):
        return x / chi0

    @staticmethod
    def f_dh(l):
        return 0.970 * (1 + (np.log(l) / 20.7)) * (1 + (np.log(l) / 22.7))

    @staticmethod
    def t(pv: float, p1v1: float, **kwargs):
        db = kwargs.get('db')
        material = kwargs.get('material')
        es = 14.1  # MeV
        chi0 = radiation_length(db=db, material=material)
        x = thickness(material, k_out=pv, k_in=p1v1, db=db)
        return DifferentialHighland.f_dh(DifferentialHighland.l(x, chi0)) * (es / pv) ** 2 * (1 / chi0)


class ICRU:
    """"""
    @staticmethod
    def t(pv: float, p1v1: float, **kwargs):
        pass


class ICRUProtons:
    """"""
    @staticmethod
    def t(pv: float, p1v1: float, **kwargs):
        db = kwargs.get('db')
        material = kwargs['material']
        es = 15.0  # MeV
        chi_s = scattering_length(material=material, db=db)
        return (es / pv) ** 2 * (1 / chi_s)


class DifferentialMoliere:
    """"""
    @staticmethod
    def t(pv: float, p1v1: float, **kwargs):
        db = kwargs.get('db')
        material = kwargs['material']
        es = 15.0  # MeV
        chi_s = scattering_length(material=material, db=db)
        return DifferentialMoliere.f_dm(p1v1, pv) * (es / pv) ** 2 * (1 / chi_s)

    @staticmethod
    def f_dm(p1v1: float, pv: float):
        if pv <= 0:
            raise ValueError("'pv' must be > 0.")
        if p1v1 <= 0:
            raise ValueError("'p1v1' must be > 0.")
        if p1v1 <= pv:
            raise ValueError("Initial 'p1v1' must be larger than final 'pv'.")
        return 0.5244 \
               + 0.1975 * np.log10(1 - (pv / p1v1) ** 2) \
               + 0.2320 * np.log10(pv) \
               - 0.0098 * np.log10(pv) * np.log10(1 - (pv / p1v1) ** 2)
