import re
from typing import Any
from pydantic import BaseModel, RootModel


def getattr_in(obj, names) -> Any:
    """
    Calls getattr on obj and the successive returns from getattr
    using names as attribute names in order. If any return is None,
    terminate early, and return None.
    Like clojure (get-in obj names) or (some-> obj .name1 .name2 ...)
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


def is_identifiable(o: Any) -> bool:
    """
    Determine if object is identifiable. An object is considered identifiable if
    contains a `ga4gh_digest` attribute

    :param o: Object
    :return: `True` if `o` has `ga4gh_digest` attribute. `False` otherwise.
    """
    return getattr_in(o, ['ga4gh', 'identifiable'])


def is_literal(o: Any) -> bool:
    return isinstance(o, (str, int, float, complex, bool))


def is_list(o: Any) -> bool:
    return isinstance(o, list)


def is_curie_type(o: Any) -> bool:
    """
    Returns true if the object is a str-like matching the CURIE pattern.
    If object is a Pydantic custom root type, extracts the value first,
    which enables (for example) a GA4GH IRI pydantic model object to be passed.
    """
    if isinstance(o, RootModel):
        o = o.root
    return re.match(r'[a-zA-Z0-9.]+:\S+', o)


def is_pydantic_instance(o: Any) -> bool:
    return isinstance(o, BaseModel)


def pydantic_copy(obj: BaseModel) -> BaseModel:
    pydantic_class = type(obj)
    if not issubclass(pydantic_class, BaseModel):
        raise RuntimeError("Argument was not a pydantic model: " + str(pydantic_class))

    # Treat RootModel differently, it's a thin wrapper of another object, has no fields
    if issubclass(pydantic_class, RootModel):
        return pydantic_class.model_construct(obj.model_dump())
    else:
        return pydantic_class.model_construct(**obj.model_dump())
