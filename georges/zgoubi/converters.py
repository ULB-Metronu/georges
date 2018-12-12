import zgoubidoo
from zgoubidoo import ureg


class Converters:
    @staticmethod
    def quadrupole(e):
        return zgoubidoo.commands.Quadrupole(e.name,
                                             XL=e.LENGTH * ureg.m,
                                             B0=0.1 * ureg.tesla,
                                             XCE=0.0 * ureg.centimeter,
                                             YCE=0.0 * ureg.centimeter,
                                             ALE=0. * ureg.degree,
                                             XPAS=1 * ureg.mm,
                                             KPOS=1,
                                             )

    @staticmethod
    def sbend(e):
        return zgoubidoo.commands.Dipole(e.name,
                                         RM=200 * ureg.cm,
                                         AT=e.ANGLE * ureg.radian,
                                         B0=10.745 * ureg.kilogauss,
                                         N=0,
                                         IL=2,
                                         ACENT=30 * ureg.degree / 2,
                                         IORDRE=25,
                                         KPOS=2,
                                         XPAS=0.1 * ureg.centimeter,
                                         OMEGA_E=30 * ureg.degree / 2,
                                         OMEGA_S=-30 * ureg.degree / 2,
                                         R1_E=1e9 * ureg.centimeter,
                                         U1_E=1e9 * ureg.centimeter,
                                         U2_E=1e9 * ureg.centimeter,
                                         R2_E=1e9 * ureg.centimeter,
                                         RE=200 * ureg.centimeter,
                                         TE=0 * ureg.radian,
                                         RS=200 * ureg.centimeter,
                                         TS=0 * ureg.radian
                                         )

    @staticmethod
    def collimator(e):
        return zgoubidoo.commands.Drift(e.name,
                                        XL=e.LENGTH * ureg.m,
                                        )
