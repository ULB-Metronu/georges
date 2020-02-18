"""
TODO
"""
from __future__ import annotations
from typing import Optional, Any, Tuple, Dict, Callable
import inspect
import uuid
import numpy as _np
from pint import UndefinedUnitError as _UndefinedUnitError
from ... import ureg as _ureg
from ..integrators import IntegratorType, MadXIntegrator
from ..apertures import circular_aperture_check, \
    rectangular_aperture_check, \
    elliptical_aperture_check, \
    phase_space_aperture_check
from georges_core.patchable import Patchable as _Patchable


class ManzoniException(Exception):
    """Exception raised for errors in the Manzoni elements module."""

    def __init__(self, m):
        self.message = m


class ManzoniAttributeException(ManzoniException):
    """Exception raised for errors in the Manzoni elements module."""
    pass


class ElementType(type):
    """
    Dark magic.
    Be careful.

    TODO
    """

    def __new__(mcs, name: str, bases: Tuple[ElementType, type, ...], dct: Dict[str, Any]):
        # Insert a default initializer (constructor) in case one is not present
        if '__init__' not in dct:
            def default_init(self,
                             label1: str = '',
                             integrator: IntegratorType = MadXIntegrator,
                             *params, **kwargs):
                """Default initializer for all Commands."""
                defaults = {}
                if 'post_init' in dct:
                    defaults = {
                        parameter_name: parameter_value.default
                        for parameter_name, parameter_value in inspect.signature(dct['post_init']).parameters.items()
                        if parameter_value.default is not inspect.Parameter.empty
                    }
                bases[0].__init__(self,
                                  label1,
                                  integrator,
                                  dct.get('PARAMETERS', {}),
                                  *params, **{**defaults, **kwargs})
                if 'post_init' in dct:
                    dct['post_init'](self, **kwargs)

            dct['__init__'] = default_init

        # Collect all post_init arguments
        if '_POST_INIT' not in dct:
            dct['_POST_INIT'] = {}
        if 'post_init' in dct and len(bases) > 0:
            dct['_POST_INIT'] = {*getattr(bases[0], '_POST_INIT', {}), *dct['post_init'].__code__.co_varnames}

        # Add PARAMETERS from the base class
        try:
            dct['PARAMETERS'] = {**getattr(bases[0], 'PARAMETERS', {}), **dct.get('PARAMETERS', {})}
        except IndexError:
            pass

        return super().__new__(mcs, name, bases, dct)

    def __init__(cls, name: str, bases: Tuple[type, ...], dct: Dict[str, Any]):
        super().__init__(name, bases, dct)
        if cls.__doc__ is not None:
            cls.__doc__ = cls.__doc__.rstrip()
            cls.__doc__ += """

    .. rubric:: Command attributes

    Attributes:
            """
            for k, v in cls.PARAMETERS.items():
                if isinstance(v, tuple) and len(v) >= 2:
                    cls.__doc__ += f"""
        {k}='{v[0]}' ({type(v[0]).__name__}): {v[1]}
            """

    def __getattr__(cls, key: str):
        try:
            if key.endswith('_'):
                return cls.PARAMETERS[key.rstrip('_')][2]
            else:
                return cls.PARAMETERS[key][0]
        except KeyError:
            raise AttributeError(key)

    def __getitem__(cls, item: str):
        try:
            return cls.PARAMETERS[item]
        except KeyError:
            raise AttributeError(item)

    def __contains__(cls, item) -> bool:
        return item in cls.PARAMETERS


class Element(metaclass=ElementType):
    """Test test test.

    More info on this wonderful class.
    TODO
    """
    PARAMETERS: dict = {
        'LABEL1': ('', 'Primary label for the Zgoubi command (default: auto-generated hash).'),
    }
    """Parameters of the element, with their default value and their description ."""

    def __init__(self, label1: str = '', *params, **kwargs):
        """
        TODO
        Args:
            label1:
            label2:
            *params:
            **kwargs:
        """
        self._attributes = {}
        for d in (Element.PARAMETERS,) + params:
            self._attributes = dict(self._attributes, **{k: v[0] for k, v in d.items()})
        for k, v in kwargs.items():
            if k not in self._POST_INIT:
                setattr(self, k, v)
        if label1:
            self._attributes['LABEL1'] = label1
        if not self._attributes['LABEL1']:
            self.generate_label()
        Element.post_init(self, **kwargs)

    def generate_label(self, prefix: str = ''):
        """

        Args:
            prefix:

        Returns:

        """
        self._attributes['LABEL1'] = '_'.join(filter(None, [
            prefix,
            str(uuid.uuid4().hex)
        ]))[:20]
        return self

    def post_init(self, **kwargs):  # -> NoReturn:
        """
        TODO
        Args:
            **kwargs: all arguments from the initializer (constructor) are passed to ``post_init`` as keyword arguments.

        """
        pass

    def __getattr__(self, a: str) -> Any:
        """
        TODO
        Args:
            a:

        Returns:

        """
        if self._attributes.get(a) is None:
            try:
                return super().__getattribute__(a)
            except AttributeError:
                return None
        return self._attributes[a]

    def __setattr__(self, k: str, v: Any):
        """
        Custom attribute setter; all non-protected (starting with a '_') attributes in upper-case are considered
        parameters of the `Command`. As such, the method will verify that they are indeed part of the command
        definition, if not an exception is raised. For valid attributes, their dimensionality is verified against the
        command definition (the dimension of the parameter's default value).

        It is also possible to use unit inference by appending an underscore to the attributes' name. In that
        case the unit of the default value is implicitely used. This is useful in case it is known that the parameter's
        numerical value is expressed in Zgoubi's default units set.

        Examples:
            >>> c = Command()
            >>> c.LABEL1 = 'FOOBAR'

        Args:
            k: a string representing the attribute
            v: the attribute's value to be set

        Raises:
            A ZgoubidooException is raised in case the parameter is not part of the class definition or if it has
            invalid dimension.
        """
        if k.startswith('_') or not k.isupper():
            super().__setattr__(k, v)
        else:
            k_ = k.rstrip('_')
            if k_ not in self._attributes.keys():
                raise ManzoniAttributeException(f"The parameter {k_} is not part of the {self.__class__.__name__} "
                                                  f"definition.")

            default = self._retrieve_default_parameter_value(k_)
            try:  # Avoid a bug in pint where a string starting with '#' cannot be parsed
                default = default.lstrip('#')
            except AttributeError:
                pass
            if isinstance(v, (int, float)) and k.endswith('_'):
                v = _ureg.Quantity(v, _ureg.Quantity(default).units)
            elif isinstance(v, str) and default is not None and not isinstance(default, str) and not v.startswith('#'):
                try:
                    v = _ureg.Quantity(v)
                except _UndefinedUnitError:
                    pass
            try:
                dimension = v.dimensionality
            except AttributeError:
                dimension = _ureg.Quantity(1).dimensionality  # No dimension
            try:
                if default is not None and dimension != _ureg.Quantity(default).dimensionality:
                    raise ManzoniAttributeException(f"Invalid dimension ({dimension} "
                                                    f"instead of {_ureg.Quantity(default).dimensionality}) "
                                                    f"for parameter {k_}={v} of {self.__class__.__name__}."
                                                    )
            except (ValueError, TypeError, _UndefinedUnitError):
                pass
            self._attributes[k_] = v

    def _retrieve_default_parameter_value(self, k: str) -> Any:
        """
        Retrieve the default value of a given parameter as defined in the Command definition (class hierarchy).

        Args:
            k: the parameter for which the default value is requested.

        Returns:
            the default value of the Command's parameter 'k'.
        """
        try:
            return self.PARAMETERS[k][0]
        except (TypeError, IndexError):
            return self.PARAMETERS[k]

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        """
        Provides the string representation of the command in the Zgoubi input file format.

        Returns:
            The string representation.

        Examples:
            >>> c = Command('my_label_1', 'my_label_2')
            >>> str(c)
        """
        return f"{self.__class__.__name__}: {str(self.attributes)}"

    @property
    def attributes(self) -> Dict[str, _ureg.Quantity]:
        """All attributes.

        Provides a dictionary with all attributes for the command.

        Returns: dictionnary with all attributes.

        """
        return self._attributes

    @property
    def defaults(self) -> Dict[str, _ureg.Quantity]:
        """Default attributes.

        Provides a dictionary with all attributes that have been assigned a default value.

        Returns: dictionary with all default attributes.
        """
        return {k: v for k, v in self._attributes.items() if v == self.PARAMETERS.get(k)[0]}

    @property
    def nondefaults(self) -> Dict[str, _ureg.Quantity]:
        """Non default attributes.

        Provides a dictionary with all attributes that have been assigned a non default value.

        Returns: dictionary with all non default attributes.
        """
        return {k: v for k, v in self._attributes.items() if v != self.PARAMETERS.get(k)[0]}


class ManzoniElement(Element, _Patchable):

    def __init__(self,
                 label1: str = '',
                 integrator: IntegratorType = MadXIntegrator,
                 *params,
                 **kwargs
                 ):
        """

        Args:
            label1:
            integrator:
            *params:
            **kwargs:
        """
        super().__init__(label1, *params, **kwargs)
        self._integrator = integrator
        self._cache: Optional[Any] = None
        self._frozen: bool = False

    def propagate(self,
                  beam: _np.ndarray,
                  out: Optional[_np.ndarray] = None,
                  global_parameters: list = None,
                  ) -> Tuple[_np.ndarray, _np.ndarray]:
        """

        Args:
            beam:
            out:
            global_parameters:

        Returns:

        """

        return self.integrator.propagate(self, beam, out, global_parameters)

    def check_aperture(self,
                       beam: _np.ndarray,
                       out: _np.ndarray
                       ):
        """

        Args:
            beam:
            out:

        Returns:

        """
        if self.APERTYPE is not None:
            aperture_check, aperture_parameters = self.aperture
            out = _np.compress(
                aperture_check(out, aperture_parameters),
                out,
                axis=0,
            )
        return beam, out

    def freeze(self):
        """

        Returns:

        """
        self.cache  # Calls it!
        self._frozen = True
        return self

    def unfreeze(self):
        """

        Returns:

        """
        self._frozen = False
        return self

    @property
    def frozen(self):
        """

        Returns:

        """
        return self._frozen

    @property
    def unfrozen(self):
        """

        Returns:

        """
        return not self._frozen

    @property
    def integrator(self) -> IntegratorType:
        """

        Returns:

        """
        return self._integrator

    @integrator.setter
    def integrator(self, integrator: IntegratorType):
        self._integrator = integrator

    @property
    def parameters(self) -> list:
        return []

    @property
    def aperture(self) -> Optional[Tuple[Callable, _np.ndarray]]:
        if self.APERTYPE.upper() == 'CIRCULAR':
            return (
                circular_aperture_check,
                _np.array([
                    self.APERTURE[0].m_as('meter')
                ])
            )
        elif self.APERTYPE.upper() == 'RECTANGULAR':
            return (
                rectangular_aperture_check,
                _np.array([
                    self.APERTURE[0].m_as('meter'), self.APERTURE[1].m_as('meter')
                ])
            )
        elif self.APERTYPE.upper() == 'ELLIPTICAL':
            return (
                elliptical_aperture_check,
                _np.array([
                    self.APERTURE[0].m_as('meter'), self.APERTURE[1].m_as('meter')
                ])
            )
        elif self.APERTYPE.upper() == 'PHASE':
            return (
                phase_space_aperture_check,
                _np.array([
                    self.APERTURE[0], self.APERTURE[1], self.APERTURE[2], self.APERTURE[3]
                ])
            )
        else:
            return None

    def clear_cache(self):
        self._cache = None

    @property
    def cache(self) -> list:
        if not self.frozen:
            self._cache = self.integrator.cache(self)
        return self._cache
