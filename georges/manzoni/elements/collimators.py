"""
TODO
"""
from ... import ureg as _ureg
from .magnets import Drift as _Drift


class Collimator(_Drift):
    """
    Define a Collimator.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Collimator with their default values.

    Examples:
        >>> c1 = Collimator('C1', L=10*_ureg.cm, APERTYPE='ELIPTICAL', APERTURE=[2*_ureg.cm, 3*_ureg.cm])
        >>> c1 #doctest: +NORMALIZE_WHITESPACE
            Collimator: {'NAME': 'C1',
                         'AT_ENTRY': <Quantity(0, 'meter')>,
                         'AT_CENTER': <Quantity(0, 'meter')>,
                         'AT_EXIT': <Quantity(0, 'meter')>,
                         'L': <Quantity(10, 'centimeter')>,
                         'APERTYPE': 'ELIPTICAL',
                         'APERTURE': [<Quantity(2, 'centimeter')>, <Quantity(3, 'centimeter')>]}
    """

    PARAMETERS = {
        "APERTYPE": (None, "Aperture type (CIRCULAR, ELLIPTICAL, RECTANGULAR or PHASE_SPACE)"),
        "APERTURE": ([None, None, None, None], ""),
    }


class CircularCollimator(Collimator):
    """
    Define a CircularCollimator.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the CircularCollimator with their default values.

    Examples:
        >>> c2 = CircularCollimator('C2', L=10*_ureg.cm, APERTURE=[2*_ureg.cm])
        >>> c2 #doctest: +NORMALIZE_WHITESPACE
            CircularCollimator: {'NAME': 'C2',
                         'AT_ENTRY': <Quantity(0, 'meter')>,
                         'AT_CENTER': <Quantity(0, 'meter')>,
                         'AT_EXIT': <Quantity(0, 'meter')>,
                         'L': <Quantity(10, 'centimeter')>,
                         'APERTYPE': 'CIRCULAR',
                         'APERTURE': [<Quantity(2, 'centimeter')>]}
    """

    PARAMETERS = {
        "APERTYPE": ("CIRCULAR", "Aperture type"),
    }


class EllipticalCollimator(Collimator):
    """
    Define an EllipticalCollimator.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the EllipticalCollimator with their default values.

    Examples:
        >>> c2 = EllipticalCollimator('C2', L=10*_ureg.cm, APERTURE=[2*_ureg.cm, 3*_ureg.cm])
        >>> c2 #doctest: +NORMALIZE_WHITESPACE
            EllipticalCollimator: {'NAME': 'C2',
                         'AT_ENTRY': <Quantity(0, 'meter')>,
                         'AT_CENTER': <Quantity(0, 'meter')>,
                         'AT_EXIT': <Quantity(0, 'meter')>,
                         'L': <Quantity(10, 'centimeter')>,
                         'APERTYPE': 'ELLIPTICAL',
                         'APERTURE': [<Quantity(2, 'centimeter')>, <Quantity(3, 'centimeter')>]}
    """

    PARAMETERS = {
        "APERTYPE": ("ELLIPTICAL", "Aperture type"),
    }


class RectangularCollimator(Collimator):
    """
    Define a RectangularCollimator.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the RectangularCollimator with their default values.

    Examples:
        >>> c3 = RectangularCollimator('C3', L=10*_ureg.cm, APERTURE=[2*_ureg.cm, 3*_ureg.cm])
        >>> c3 #doctest: +NORMALIZE_WHITESPACE
            RectangularCollimator: {'NAME': 'C3',
                         'AT_ENTRY': <Quantity(0, 'meter')>,
                         'AT_CENTER': <Quantity(0, 'meter')>,
                         'AT_EXIT': <Quantity(0, 'meter')>,
                         'L': <Quantity(10, 'centimeter')>,
                         'APERTYPE': 'RECTANGULAR',
                         'APERTURE': [<Quantity(2, 'centimeter')>, <Quantity(3, 'centimeter')>]}
    """

    PARAMETERS = {
        "APERTYPE": ("RECTANGULAR", "Aperture type"),
    }


class Dump(RectangularCollimator):
    """
    Define a Dump.

    Attributes:
        PARAMETERS (dict): Dictionary containing the parameters of the Dump with their default values.

    Examples:
        >>> d1 = Dump('D1', L=10*_ureg.cm)
        >>> d1 #doctest: +NORMALIZE_WHITESPACE
            Dump: {'NAME': 'D1',
                   'AT_ENTRY': <Quantity(0, 'meter')>,
                   'AT_CENTER': <Quantity(0, 'meter')>,
                   'AT_EXIT': <Quantity(0, 'meter')>,
                   'L': <Quantity(10, 'centimeter')>,
                   'APERTYPE': 'RECTANGULAR',
                   'APERTURE': <Quantity(0.0, 'centimeter')>}
    """

    PARAMETERS = {
        "APERTURE": [0.0 * _ureg.cm, 0.0 * _ureg.cm],
    }
