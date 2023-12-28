from importlib.metadata import version

__version__ = version("georges")

from georges_core import Q_, Kinematics, particles, ureg
from georges_core.distribution import *
from georges_core.madx import madx
from georges_core.sequences import *
from georges_core.sequences.betablock import BetaBlock

from .vis import ManzoniMatplotlibArtist, ManzoniPlotlyArtist
