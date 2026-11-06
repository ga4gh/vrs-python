"""decorators for vrs-python"""

from collections.abc import Callable
from typing import TypeVar

_Instance = TypeVar("_Instance")
_Value = TypeVar("_Value")


def lazy_property(fn: Callable[[_Instance], _Value]) -> property:
    """Provide a decorator that makes a property lazy-evaluated.

    :param fn: Function whose result is cached on each instance.
    :return: Property that caches its value on first access.
    """
    attr_name = "_lazy_" + fn.__name__

    @property
    def _lazy_property(self: _Instance) -> _Value:
        """Compute and cache the property value on first access.

        :return: Cached or newly computed value.
        """
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazy_property
