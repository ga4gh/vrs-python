import re
from typing import Any
from pydantic import BaseModel


def is_identifiable(o: Any) -> bool:
    """
    Determine if object is identifiable. An object is considered identifiable if
    contains a `ga4gh_digest` attribute

    :param o: Object
    :return: `True` if `o` has `ga4gh_digest` attribute. `False` otherwise.
    """
    return hasattr(o, "ga4gh_digest")


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
    if isinstance(o, BaseModel) and hasattr(o, "__root__"):
        o = o.__root__
    return re.match(r'[a-zA-Z0-9.]+:\S+', o)




def is_pydantic_instance(o: Any) -> bool:
    return isinstance(o, BaseModel)


def pydantic_copy(obj: BaseModel) -> BaseModel:
    pydantic_class = type(obj)
    if not issubclass(pydantic_class, BaseModel):
        raise RuntimeError("Argument was not a pydantic model: " + str(pydantic_class))
    return pydantic_class(**obj.dict())
