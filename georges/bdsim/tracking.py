from .. import beamline
from .bdsim import BDSim
from .. import physics


class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def track(**kwargs):
    """Compute the distribution of the beam as it propagates through the beamline..
    :param kwargs: parameters are:
        - line: the beamline on which twiss will be run
        - context: the associated context on which MAD-X is run
    """
    # Process arguments
    line = kwargs.get('line', None)
    b = kwargs.get('beam', None)
    context = kwargs.get('context', None)
    options = kwargs.get('options', None)

    if line is None or b is None or context is None:
        raise TrackException("Beamline, Beam, Context and BDsim objects need to be defined.")

    if options is None:
        print("Warning : no options is provided")

    bd = BDSim(beamlines=[line], **kwargs)

    # Write the input file for bdsim
    bd_beam = b.distribution.copy()
    p0 = physics.energy_to_momentum(b.energy)

    bd.track(bd_beam, p0, **kwargs)

    # Add options for Bdsim
    bd.set_options(options)

    # Write bdsim input
    errors = bd.run(**kwargs).fatals

