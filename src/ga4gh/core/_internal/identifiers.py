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
from typing import Union
from pydantic import BaseModel
from canonicaljson import encode_canonical_json

from .digests import sha512t24u
from .pydantic import (
    is_pydantic_instance,
    is_curie_type,
    is_identifiable,
    getattr_in,
    get_pydantic_root,
    is_pydantic_custom_str_type
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
    return str(ir).startswith(ns_w_sep)


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


def ga4gh_identify(vro, type_prefix_map=None):
    """return the GA4GH digest-based id for the object, as a CURIE
    (string).  Returns None if object is not identifiable.

    >>> import ga4gh.vrs
    >>> ival = ga4gh.vrs.models.SimpleInterval(start=44908821, end=44908822)
    >>> location = ga4gh.vrs.models.Location(sequence_id="ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", interval=ival)
    >>> ga4gh_identify(location)
    'ga4gh:VSL.u5fspwVbQ79QkX6GHLF8tXPCAXFJqRPx'

    """
    digest = ga4gh_digest(vro)
    pfx = vro.prefix
    ir = f"{namespace}{curie_sep}{pfx}{ref_sep}{digest}"
    return ir


def ga4gh_digest(vro, do_compact=True):
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
    print("ga4gh_digest: vro: " + str(vro))
    s = ga4gh_serialize(vro, do_compact=do_compact)
    return sha512t24u(s)


def export_pydantic_model(obj, exclude_none=True):
    # Export Pydantic model to raw Python object. If a custom root object,
    # return that. If a Model with fields, exports to a dict.
    # TODO maybe just call export_pydantic_model on the .__root__ instead of
    # get_pydantic_root, taking whatever it is. Recursion should terminate fine.
    if isinstance(obj, BaseModel):
        # try custom root type first, if not, assume it's a normal class
        obj2 = get_pydantic_root(obj)
        if obj2 != obj:
            obj = obj2
        else:
            obj = obj.model_dump(exclude_none=exclude_none)
    return obj


"""
TODO: discussed making all objects digestible. If no digest keys defined,
include all fields. We first need to define keys for all model objects.
"""


def _serialize_compact(
    input_obj: Union[BaseModel, dict, str],
    enref=True
) -> Union[str, dict]:
    """
    Compacts a Pydantic model prior to serialization as canonicaljson.
    Rolls up nested identifiable objeccts and shortens
    ga4gh identifier custom types to their digest segments.
    """
    if input_obj is None:
        return None
    output_obj = input_obj
    if is_pydantic_custom_str_type(input_obj):
        val = export_pydantic_model(input_obj)
        if is_curie_type(val) and is_ga4gh_identifier(val):
            val = val[len("ga4gh:"):]
            if "." in val:
                val = val[val.index(".") + 1:]
        output_obj = val
    elif is_pydantic_instance(input_obj):
        # Take static key set from the object, or use all fields
        include_keys = getattr_in(input_obj, ["ga4gh", "keys"])
        if include_keys is None or len(include_keys) == 0:
            exported = export_pydantic_model(input_obj)
            include_keys = exported.keys()
        # Serialize each field value
        output_obj = {
            k: _serialize_compact(getattr(input_obj, k), enref=True)
            for k in include_keys
            if hasattr(input_obj, k)  # check if None?
        }
        # Collapse identifiable object into digest
        # TODO store the object itself as well, with the digest,
        # so caller has this information and doesn't need to recompute.
        if is_identifiable(input_obj) and enref:
            # Collapse the obj into its digest, not re-running compaction
            output_obj = ga4gh_digest(output_obj, do_compact=False)
    else:
        exported_obj = export_pydantic_model(input_obj)
        if type(exported_obj) in [list, set]:
            output_obj = [_serialize_compact(elem) for elem in exported_obj]
    return output_obj


def identify_all(
    input_obj: Union[BaseModel, dict, str]
) -> Union[str, dict]:
    """
    Compacts a Pydantic model prior to serialization as canonicaljson.
    Rolls up nested identifiable objeccts and shortens
    ga4gh identifier custom types to their digest segments.
    """
    if input_obj is None:
        return None
    output_obj = input_obj
    print("input_obj: " + str(input_obj))
    if is_pydantic_custom_str_type(input_obj):
        val = export_pydantic_model(input_obj)
        if is_curie_type(val) and is_ga4gh_identifier(val):
            val = parse_ga4gh_identifier(val)["digest"]
        output_obj = val
    elif is_pydantic_instance(input_obj):
        # Take static key set from the object, or use all fields
        include_keys = getattr_in(input_obj, ["ga4gh", "keys"])
        # TODO Add keys to each Model class
        if include_keys is None or len(include_keys) == 0:
            include_keys = export_pydantic_model(input_obj).keys()
        # Serialize each field value
        output_obj = {
            k: identify_all(getattr(input_obj, k))
            for k in include_keys
            if hasattr(input_obj, k)  # check if None?
        }
        print("output_obj: " + str(output_obj))
        # Collapse identifiable object into digest
        # TODO store the object itself as well, with the digest,
        # so caller has this information and doesn't need to recompute.
        if is_identifiable(input_obj):
            # Compute digest for updated object, not re-running compaction
            output_obj["digest"] = ga4gh_digest(output_obj, do_compact=False)
    else:
        exported_obj = export_pydantic_model(input_obj)
        if type(exported_obj) in [list, set]:
            output_obj = [_serialize_compact(elem) for elem in exported_obj]
    return output_obj


def ga4gh_serialize(obj, do_compact=True):
    if do_compact:
        obj = _serialize_compact(obj, enref=False)
    if isinstance(obj, dict):
        return encode_canonical_json(obj)
    else:
        return obj.encode("utf-8")


def scrape_model_metadata(obj, meta={}) -> dict:
    """
    """
    assert isinstance(obj, BaseModel)
    # map of class names to 'identifiable' and 'keys'
    meta = {}
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
        # RE
    return meta
