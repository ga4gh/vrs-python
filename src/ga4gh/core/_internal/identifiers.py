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

import re
from typing import Union
from pydantic import BaseModel
from canonicaljson import encode_canonical_json

from .digests import sha512t24u
from .pydantic import (
    is_list,
    is_pydantic_instance,
    is_curie_type,
    is_identifiable,
    is_literal,
    getattr_in)

__all__ = "ga4gh_digest ga4gh_identify ga4gh_serialize is_ga4gh_identifier parse_ga4gh_identifier".split()

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


def ga4gh_digest(vro):
    """return the GA4GH digest for the object

    >>> import ga4gh.vrs
    >>> ival = ga4gh.vrs.models.SimpleInterval(start=44908821, end=44908822)
    >>> location = ga4gh.vrs.models.Location(sequence_id="ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", interval=ival)
    >>> ga4gh_digest(location)
    'u5fspwVbQ79QkX6GHLF8tXPCAXFJqRPx'

    """

    assert is_identifiable(vro), "ga4gh_digest called with non-identifiable object"
    return sha512t24u(ga4gh_serialize(vro))


def ga4gh_digest2(vro, do_compact=True):
    """
    Return the GA4GH digest for the object

    >>> import ga4gh.vrs
    >>> ival = ga4gh.vrs.models.SimpleInterval(start=44908821, end=44908822)
    >>> location = ga4gh.vrs.models.Location(sequence_id="ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", interval=ival)
    >>> ga4gh_digest(location)
    'u5fspwVbQ79QkX6GHLF8tXPCAXFJqRPx'

    """

    # assert is_identifiable(vro), "ga4gh_digest called with non-identifiable object"
    s = ga4gh_serialize2(vro, do_compact=do_compact)
    return sha512t24u(s)


def export_pydantic_model(obj):
    # Export Pydantic model to raw Python object. Right now supporting
    # dict, or a custom root type that is a str
    if isinstance(obj, BaseModel):
        # try custom root type first, if not, assume it's a normal class
        obj2 = get_root_str(obj)
        if obj2 != obj:
            obj = obj2
        else:
            obj = obj.dict()
    return obj


def get_root_str(obj: Union[str, BaseModel]) -> str:
    """
    If o is a Pydantic custom root type that is a string, return that string
    """
    if (isinstance(obj, BaseModel)
            and hasattr(obj, "__root__")
            and isinstance(obj.__root__, str)):
        return obj.__root__
    return obj


def is_pydantic_custom_str_type(obj: BaseModel):
    return hasattr(obj, "__root__") and isinstance(obj.__root__, str)


def _serialize_compact(
    input_obj: Union[BaseModel, dict, str]
) -> Union[str, dict]:
    """
    Compacts a Pydantic model prior to serialization as canonicaljson.
    Rolls up nested identifiable objeccts and shortens
    ga4gh identifier custom types to their digest segments.
    """
    if input_obj is None:
        return None
    print(input_obj)
    output_obj = input_obj
    # obj = export_pydantic_model(input_obj)
    if is_pydantic_custom_str_type(input_obj):
        val = export_pydantic_model(input_obj)
        if is_curie_type(val) and is_ga4gh_identifier(val):
            val = val[len("ga4gh:"):]
            if "." in val:
                val = val[val.index(".") + 1:]
        output_obj = val
    elif is_pydantic_instance(input_obj):
        exported_obj = export_pydantic_model(input_obj)
        include_keys = getattr_in(input_obj, ["ga4gh", "keys"])
        if include_keys is None or len(include_keys) == 0:
            include_keys = exported_obj.keys()
        output_obj = {
            k: _serialize_compact(exported_obj[k])
            for k in include_keys
            if exported_obj[k] is not None
        }
        if is_identifiable(input_obj):
            # Collapse the obj into its digest, not re-running compaction
            output_obj = ga4gh_digest2(output_obj, do_compact=False)
            # output_obj = sha512t24u(encode_canonical_json(output_obj))
    else:
        exported_obj = export_pydantic_model(input_obj)
        if type(exported_obj) in [list, set]:
            output_obj = [_serialize_compact(elem) for elem in exported_obj]
    return output_obj


def ga4gh_serialize2(obj, do_compact=True):
    if do_compact:
        obj = _serialize_compact(obj)
    return encode_canonical_json(obj)


def ga4gh_serialize(vro):
    """
    Serialize object into a canonical format

    Briefly:
    * format is json
    * keys sorted in unicode order (=ascii order for our use)
    * no "insignificant" whitespace, as defined in rfc7159§2
    * MUST use two-char escapes when available, as defined in rfc7159§7
    * UTF-8 encoded
    * nested identifiable objects are replaced by their identifiers
    * arrays of identifiers are sorted lexographically

    These requirements are a distillation of several proposals which
    have not yet been ratified.

    >>> import ga4gh.vrs
    >>> ival = ga4gh.vrs.models.SimpleInterval(start=44908821, end=44908822)
    >>> location = ga4gh.vrs.models.Location(sequence_id="ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", interval=ival)
    >>> ga4gh_serialize(location)
    b'{"interval":{"end":44908822,...,"type":"SequenceLocation"}'

    """

    def dictify(vro, enref=True):
        """
        Recursively converts (any) object to dictionary prior to
        serialization

        enref: if True, replace nested identifiable objects with
        digests ("enref" is opposite of "de-ref")
        """

        if vro is None:    # pragma: no cover
            return None

        if is_literal(vro):
            v = vro._value
            if is_curie_type(vro):
                if is_ga4gh_identifier(v):
                    # CURIEs are stripped to just the digest so that digests are independent of type prefixes
                    v = v.split(ref_sep, 1)[1]
            return v

        if isinstance(vro, str):
            v = vro
            if is_ga4gh_identifier(v):
                v = v.split(ref_sep, 1)[1]
            return v

        if is_pydantic_instance(vro):
            if is_identifiable(vro) and enref:
                return ga4gh_digest(vro)

            d = {k: dictify(vro[k], enref=True)
                 for k in vro
                 if k in vro.ga4gh_digest.keys}
            return d

        if is_list(vro):
            if is_curie_type(vro[0]):
                return sorted(dictify(o) for o in vro.data)
            return sorted([dictify(o) for o in vro.typed_elems])

        raise ValueError(f"Don't know how to serialize {vro}")    # pragma: no cover

    # The canonicaljson package does everything we want. Use that with
    # the hope that it will be upward compatible with a future
    # ratified proposal for json canonicalization.
    #
    # The following alternative does the same thing for our use case.
    # It's included here as an outline for anyone implementing in
    # another language.  (canonicaljson escapes unicode characters, as
    # required by VRS, but this doesn't apply to any known uses so
    # these are equivalent.)

    # >> import json
    # >> def cjdump(a):
    # >>     return json.dumps(a, sort_keys=True, separators=(',',':'),
    #                          indent=None).encode("utf-8")

    vro_dict = dictify(vro, enref=False)
    return encode_canonical_json(vro_dict)
