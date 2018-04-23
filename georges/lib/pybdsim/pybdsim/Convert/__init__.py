"""
Module for various conversions.

"""

from ._MadxTfs2Gmad import MadxTfs2Gmad
from ._MadxTfs2GmadStrength import MadxTfs2GmadStrength

try:
    from _Mad8Saveline2Gmad import Mad8Saveline2Gmad
    from _Mad8Twiss2Gmad import Mad8Twiss2Gmad
    from _Mad8Twiss2Gmad import Mad8MakeOptions
    from _Mad8Twiss2Gmad import Mad8MakeApertureTemplate
    from _Mad8Twiss2Gmad import Mad8MakeCollimatorTemplate
except ImportError:
    import warnings
    msg = "Missing pymad8 dependency.  MAD8 conversion facilities excluded."
    warnings.warn(msg)
    del warnings

try:
    import pysad
    from _SadFlat2Gmad import SadFlat2GMad
except ImportError:
    pass
    #import warnings
    #msg = "Missing pysad dependency. SAD conversion facilities excluded."
    #warnings.warn(msg)
    #del warnings

try:
    from _Transport2Gmad import Transport2Gmad
except ImportError:
    import warnings
    msg = "Missing pytransport dependency.  TRANSPORT conversion facilities excluded."
    warnings.warn(msg)
    del warnings

from ._BdsimPrimaries2Inrays import BdsimPrimaries2Ptc
from ._BdsimPrimaries2Inrays import BdsimPrimaries2Madx
from ._BdsimPrimaries2Inrays import BdsimPrimaries2Mad8
