"""Provide helpful methods related to pydantic objects."""

import re
from collections.abc import Sequence
from typing import Literal, TypeVar

from pydantic import BaseModel, RootModel

_Model = TypeVar("_Model", bound=BaseModel)


def getattr_in(obj: object, names: Sequence[str]) -> object:
    """Traverse successive attributes by name.

    Stop and return None when an attribute is missing or its value is None.

    Like clojure (get-in obj names) or (some-> obj .name1 .name2 ...)

    :param obj: Object to inspect.
    :param names: Attribute names to traverse in order.
    :return: Final attribute value, or None if traversal stops early or names is empty.
    """
    names_i = 0
    v = None
    while names_i < len(names):
        v = getattr(obj, names[names_i], None)
        if v is None:
            break
        names_i += 1
        obj = v
    return v


def is_curie_type(o: object) -> re.Match[str] | Literal[False] | None:
    """Match a string or root-model value against the CURIE pattern.

    If object is a Pydantic custom root type, extracts the value first, which enables
    (for example) a GA4GH iriReference pydantic model object to be passed.

    :param o: Object to process.
    :return: Regex match for a CURIE, None for a nonmatching string, or False for other values.
    """
    if isinstance(o, RootModel):
        o = o.root
    if isinstance(o, str):
        return re.match(r"[a-zA-Z0-9.]+:\S+", o)
    return False


def is_pydantic_instance(o: object) -> bool:
    """Check whether a value is a Pydantic model.

    :param o: Object to process.
    :return: Whether the value is a Pydantic model.
    """
    return isinstance(o, BaseModel)


def get_pydantic_root(obj: object) -> object:
    """If o is a Pydantic custom root type, return the root object, else return the input obj

    :param obj: Object to inspect.
    :return: Root value for a root model, otherwise the input.
    """
    if isinstance(obj, RootModel):
        return obj.root
    return obj


def pydantic_copy(obj: _Model) -> _Model:
    """Reconstruct a model from its serialized field values.

    :param obj: Object to inspect.
    :return: New instance reconstructed from the model data.
    :raises RuntimeError: If the input is not a Pydantic model.
    """
    pydantic_class = type(obj)
    if not issubclass(pydantic_class, BaseModel):
        msg = f"Argument was not a pydantic model: {pydantic_class!r}"
        raise RuntimeError(msg)  # noqa: TRY004

    # Treat RootModel differently, it's a thin wrapper of another object, has no fields
    if issubclass(pydantic_class, RootModel):
        return pydantic_class(obj.model_dump())
    return pydantic_class(**obj.model_dump())
