import logging

try:
    import cpymad
except ImportError:
    logging.getLogger(__name__).warning("cpymad is not present, this module is loaded with reduced functionalities.")
from .madx import MadX
