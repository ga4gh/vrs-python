"""serializes, digests, and identifies GA4GH objects

In GA4GH schemas with nested objects, serialize, digest, and identify
are entangled.

For example, here is a call path for ga4gh_identify called on an Allele:
    ga4gh_identify(allele)
    + ga4gh_digest(allele)
    ++ ga4gh_serialize(allele)
    +++ ga4gh_digest(allele.location)
    ++++ ga4gh_serialize(allele.location)
    +++ ga4gh_serialize(allele.state)

For that reason, they are implemented here in one file.

"""

import contextvars
import logging
import re
from contextlib import ContextDecorator
from enum import IntEnum
from typing import Union, Optional
from pydantic import BaseModel, RootModel

from .pydantic import get_pydantic_root

__all__ = "ga4gh_digest ga4gh_identify ga4gh_serialize is_ga4gh_identifier parse_ga4gh_identifier".split()

_logger = logging.getLogger(__name__)

CURIE_NAMESPACE = "ga4gh"
CURIE_SEP = ":"
GA4GH_PREFIX_SEP = "."

GA4GH_IR_REGEXP = re.compile(r"^ga4gh:(?P<type>[^.]+)\.(?P<digest>[0-9A-Za-z_\-]{32})$")
GA4GH_DIGEST_REGEXP = re.compile(r"^[0-9A-Za-z_\-]{32}$")

ns_w_sep = CURIE_NAMESPACE + CURIE_SEP


class VrsObjectIdentifierIs(IntEnum):
    """
    Defines the state for when the `ga4gh_identify` method should compute
    an identifier ('id' attribute) for the specified object.  The options are:
      ANY - Always compute the identifier (this is the default behavior)
      GA4GH_INVALID - Compute the identifier if it is missing or is present but syntactically invalid
      MISSING - Only compute the identifier if missing

    The default behavior is safe and ensures that the identifiers are correct, 
    but at a performance cost. Where the source of inputs to `ga4gh_identify` 
    are well controlled, for example when annotating a VCF file with VRS IDs, 
    using `MISSING` can improve performance.
    """

    ANY = 0
    GA4GH_INVALID = 1
    MISSING = 2


ga4gh_compute_identifier_when = contextvars.ContextVar("ga4gh_compute_identifier_when")


class use_ga4gh_compute_identifier_when(ContextDecorator):
    """
    Context manager that defines when to compute identifiers
    for all operations within the context.  For example:

    with use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.GA4GH_INVALID):
        VCFAnnotator(...).annotate(...)

    Or:

    @use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.GA4GH_INVALID)
    def my_method():
    """

    def __init__(self, when: VrsObjectIdentifierIs):
        self.when = when
        self.token = None

    def __enter__(self):
        self.token = ga4gh_compute_identifier_when.set(self.when)

    def __exit__(self, exc_type, exc, exc_tb):
        ga4gh_compute_identifier_when.reset(self.token)


def is_ga4gh_identifier(ir):
    """

    >>> is_ga4gh_identifier("ga4gh:SQ.0123abcd")
    True

    >>> is_ga4gh_identifier("refseq:NM_01234.5")
    False

    >>> is_ga4gh_identifier(None)
    False

    """
    return str(get_pydantic_root(ir)).startswith(ns_w_sep)


def parse_ga4gh_identifier(ir):
    """
    Parses a GA4GH identifier, returning a dict with type and digest components

    >>> parse_ga4gh_identifier("ga4gh:SQ.0123abcd")
    {'type': 'SQ', 'digest': '0123abcd'}

    >>> parse_ga4gh_identifier("notga4gh:SQ.0123abcd")
    Traceback (most recent call last):
    ...
    ValueError: notga4gh:SQ.0123abcd

    """

    try:
        return GA4GH_IR_REGEXP.match(str(ir)).groupdict()
    except AttributeError as e:
        raise ValueError(ir) from e


def ga4gh_identify(vro, in_place='default'):
    """
    Return the GA4GH digest-based id for the object, as a CURIE
    (string).  Returns None if object is not identifiable.

    This function has three options for in_place editing of vro.id:
    - 'default': the standard identifier update behavior for GA4GH
        identifiable objects, this mode will update the vro.id
        field if the field is empty
    - 'always': this will update the vro.id field any time the
        identifier is computed (compute behavior is controlled by the
        use_ga4gh_compute_identifier_when context)
    - 'never': the vro.id field will not be edited in-place,
        even when empty

    TODO update example for VRS 2.0
    >>> import ga4gh.vrs
    >>> ival = ga4gh.vrs.models.SimpleInterval(start=44908821, end=44908822)
    >>> location = ga4gh.vrs.models.Location(sequence_id="ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", interval=ival)
    >>> ga4gh_identify(location)
    'ga4gh:VSL.u5fspwVbQ79QkX6GHLF8tXPCAXFJqRPx'

    """
    if vro.is_ga4gh_identifiable():
        when_rule = ga4gh_compute_identifier_when.get(VrsObjectIdentifierIs.ANY)
        obj_id = None
        if when_rule == VrsObjectIdentifierIs.ANY:
            do_compute = True
        else:
            obj_id = getattr(vro, "id", None)
            if when_rule == VrsObjectIdentifierIs.MISSING:
                do_compute = obj_id is None or obj_id == ""
            else:  # GA4GHComputeIdentifierIs.GA4GH_INVALID
                do_compute = not vro.has_valid_ga4gh_id()

        if do_compute:
            obj_id = vro.get_or_create_ga4gh_identifier(in_place)

        return obj_id

    return None


def ga4gh_digest(vro: BaseModel, overwrite=False):
    """
    Return the GA4GH digest for the object.

    do_compact: bool - true if object compaction should be performed during serialization

    TODO update example

    >>> import ga4gh.vrs
    >>> ival = ga4gh.vrs.models.SimpleInterval(start=44908821, end=44908822)
    >>> location = ga4gh.vrs.models.Location(sequence_id="ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", interval=ival)
    >>> ga4gh_digest(location)
    'u5fspwVbQ79QkX6GHLF8tXPCAXFJqRPx'

    """
    if vro.is_ga4gh_identifiable():  # Only GA4GH identifiable objects are GA4GH digestible
        return vro.get_or_create_digest(overwrite)
    else:
        return None


def replace_with_digest(val: dict) -> Union[str, dict]:
    """
    If val has a digest computed, return it, else return val
    """
    if isinstance(val, dict) and val.get("digest", None) is not None:
        return val["digest"]
    return val


def collapse_identifiable_values(obj: dict) -> dict:
    """
    Replaces dict values with their digests if they are defined.
    Does not collapse the top level object, only objects it contains.
    """
    if isinstance(obj, dict):
        obj = {
            k: replace_with_digest(collapse_identifiable_values(obj[k]))
            for k in obj.keys()
        }
    elif isinstance(obj, list) or isinstance(obj, set):
        obj = [replace_with_digest(collapse_identifiable_values(elem)) for elem in obj]
    return obj


def ga4gh_serialize(obj: BaseModel) -> Optional[bytes]:
    """
    TODO find a way to output identify_all without the 'digest' fields on subobjects,
    without traversing the whole tree again in collapse_identifiable_values.
    """
    return obj.model_dump_json().encode("utf-8")


def export_pydantic_model(obj, exclude_none=True):
    # Export Pydantic model to raw Python object. If a custom root object,
    # return that. If a Model with fields, exports to a dict.
    # TODO maybe just call export_pydantic_model on the .root instead of
    # get_pydantic_root, taking whatever it is. Recursion should terminate fine.
    if isinstance(obj, BaseModel):
        # try custom root type first, if not, assume it's a normal class
        if isinstance(obj, RootModel):
            obj = get_pydantic_root(obj)
        else:
            obj = obj.model_dump(exclude_none=exclude_none)
    return obj


# def scrape_model_metadata(obj, meta={}) -> dict:
#     """
#     For a Pydantic object obj, pull out .ga4gh.identifiable
#     and .ga4gh.keys and put them in meta keyed by the class name of obj
#     """
#     assert isinstance(obj, BaseModel)
#     name = type(obj).__name__
#     if is_pydantic_custom_str_type(obj):
#         meta[name] = {"identifiable": False, "keys": None}
#     else:
#         meta[name] = {}
#         identifiable = getattr_in(obj, ["ga4gh", "identifiable"])
#         if identifiable:
#             meta[name]["identifiable"] = identifiable
#         keys = getattr_in(obj, ["ga4gh", "keys"])
#         if keys and len(keys) > 0:
#             meta[name]["keys"] = keys
#         # TODO recurse into fields
#     return meta
