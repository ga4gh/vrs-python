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

import logging
import re
from typing import Union, Tuple
from pydantic import BaseModel, RootModel
from canonicaljson import encode_canonical_json

from .digests import sha512t24u
from .pydantic import (
    is_pydantic_instance,
    is_curie_type,
    is_identifiable,
    getattr_in,
    get_pydantic_root,
    is_pydantic_custom_type
)

__all__ = "ga4gh_digest ga4gh_identify ga4gh_serialize is_ga4gh_identifier parse_ga4gh_identifier".split()

_logger = logging.getLogger(__name__)

namespace = "ga4gh"
curie_sep = ":"
ref_sep = "."

ga4gh_ir_regexp = re.compile(r"^ga4gh:(?P<type>[^.]+)\.(?P<digest>.+)$")

ns_w_sep = namespace + curie_sep


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
        return ga4gh_ir_regexp.match(str(ir)).groupdict()
    except AttributeError as e:
        raise ValueError(ir) from e


def ga4gh_identify(vro):
    """
    Return the GA4GH digest-based id for the object, as a CURIE
    (string).  Returns None if object is not identifiable.

    TODO update example for VRS 2.0
    >>> import ga4gh.vrs
    >>> ival = ga4gh.vrs.models.SimpleInterval(start=44908821, end=44908822)
    >>> location = ga4gh.vrs.models.Location(sequence_id="ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", interval=ival)
    >>> ga4gh_identify(location)
    'ga4gh:VSL.u5fspwVbQ79QkX6GHLF8tXPCAXFJqRPx'

    """
    if is_identifiable(vro):
        digest = ga4gh_digest(vro)
        pfx = vro.ga4gh.prefix
        ir = f"{namespace}{curie_sep}{pfx}{ref_sep}{digest}"
        return ir
    return None


def ga4gh_digest(vro: BaseModel, do_compact=True):
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
    s = ga4gh_serialize(vro)
    return sha512t24u(s)


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


def ga4gh_serialize(obj: BaseModel) -> bytes:
    """
    TODO find a way to output identify_all without the 'digest' fields on subobjects,
    without traversing the whole tree again in collapse_identifiable_values.
    """
    identified = identify_all(obj)
    if isinstance(identified, dict):
        # Replace identifiable subobjects with their digests
        collapsed = collapse_identifiable_values(identified)
        if "digest" in collapsed:
            del collapsed["digest"]
        return encode_canonical_json(collapsed)
    else:
        return identified.encode("utf-8")


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


"""
TODO: discussed making all objects digestible. If no digest keys defined,
include all fields. We first need to define keys for all model objects.
"""


def identify_all(
    input_obj: Union[BaseModel, dict, str]
) -> Union[str, dict]:
    """
    Adds digests to an identifiable Pydantic object and any identifiable Pydantic
    objects in its fields, at any depth.

    Returns the identified object tree, and the tree with identified objects
    replaced with their digests.

    TODO It would be nice to have a pydantic-agnostic version of this that just takes
    a dict for input_object, and another dict that has the identifiable+keys metadata.
    Something like scrape_model_metadata can be used to generate that metadata.
    """
    if input_obj is None:
        return None
    output_obj = input_obj
    if is_pydantic_custom_type(input_obj):
        val = export_pydantic_model(input_obj)
        if isinstance(val, str) and is_curie_type(val) and is_ga4gh_identifier(val):
            val = parse_ga4gh_identifier(val)["digest"]
        output_obj = val
    elif is_pydantic_instance(input_obj):
        exported_obj = export_pydantic_model(input_obj)
        if "digest" in exported_obj and exported_obj["digest"] is not None:
            output_obj = exported_obj
        else:
            # Take static key set from the object, or use all fields
            include_keys = getattr_in(input_obj, ["ga4gh", "keys"])
            # TODO Add keys to each Model class
            if include_keys is None or len(include_keys) == 0:
                include_keys = exported_obj.keys()
            if "digest" in include_keys:
                include_keys.remove("digest")
            # Serialize each field value
            output_obj = {
                k: identify_all(getattr(input_obj, k))
                for k in include_keys
                if hasattr(input_obj, k)  # check if None?
            }
            # Assumes any obj with 'digest' should be collapsed.
            collapsed_output_obj = collapse_identifiable_values(output_obj)
            # Add a digest to the output if it is identifiable
            if is_identifiable(input_obj):
                # Compute digest for updated object, not re-running compaction
                output_obj["digest"] = ga4gh_digest(collapsed_output_obj, do_compact=False)
    else:
        exported_obj = export_pydantic_model(input_obj)
        if type(exported_obj) in [list, set]:
            output_obj = [identify_all(elem) for elem in exported_obj]
    return output_obj


def scrape_model_metadata(obj, meta={}) -> dict:
    """
    For a Pydantic object obj, pull out .ga4gh.identifiable
    and .ga4gh.keys and put them in meta keyed by the class name of obj
    """
    assert isinstance(obj, BaseModel)
    name = type(obj).__name__
    if is_pydantic_custom_str_type(obj):
        meta[name] = {"identifiable": False, "keys": None}
    else:
        meta[name] = {}
        identifiable = getattr_in(obj, ["ga4gh", "identifiable"])
        if identifiable:
            meta[name]["identifiable"] = identifiable
        keys = getattr_in(obj, ["ga4gh", "keys"])
        if keys and len(keys) > 0:
            meta[name]["keys"] = keys
        # TODO recurse into fields
    return meta
