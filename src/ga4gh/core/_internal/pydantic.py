from pydantic import BaseModel
import re
from ga4gh.vrs.models import IRI


def is_identifiable(o: any) -> bool:
    """
    Determine if object is identifiable. An object is considered identifiable if
    contains a `ga4gh_digest` attribute

    :param o: Object
    :return: `True` if `o` has `ga4gh_digest` attribute. `False` otherwise.
    """
    return hasattr(o, "ga4gh_digest")


def is_literal(o: any) -> bool:
    return isinstance(o, (str, int, float, complex, bool))


def is_list(o: any) -> bool:
    return isinstance(o, list)


def is_curie_type(o: any) -> bool:
    # if isinstance(o, pydantic.BaseModel):
    # return isinstance(o, IRI) or re.match(r'[a-zA-Z0-9.]+:\S+')
    raise RuntimeError("Not implemented: is_curie_type")


def is_pydantic_instance(o: any) -> bool:
    return isinstance(o, BaseModel)


def pydantic_copy(obj: BaseModel) -> BaseModel:
    pydantic_class = type(obj)
    if not issubclass(pydantic_class, BaseModel):
        raise RuntimeError("Argument was not a pydantic model: " + str(pydantic_class))
    return pydantic_class(**obj.dict())
