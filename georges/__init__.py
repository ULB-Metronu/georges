__version__ = "2020.1"

try:
    from georges_core import ureg, Q_
except ModuleNotFoundError:
    # TODO error handling
    # Manipulation of physical quantities (with units, etc.)
    # https://pint.readthedocs.io/en/latest/
    from pint import UnitRegistry
    ureg = UnitRegistry()
    _Q = ureg.Quantity
    ureg.define('electronvolt = e * volt = eV')
    ureg.define('electronvolt_per_c = eV / c = eV_c')
    ureg.define('electronvolt_per_c2 = eV / c**2 = eV_c2')

try:
    from georges_core import Kinematics
    from georges_core.sequences import *
    from georges_core.distribution import *
    from georges_core import particles
    from georges_core.sequences.betablock import BetaBlock
    from georges_core.madx import madx
except ModuleNotFoundError:
    pass
