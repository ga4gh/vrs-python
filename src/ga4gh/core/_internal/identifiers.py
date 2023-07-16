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

from canonicaljson import encode_canonical_json

from .digests import sha512t24u
from .pydantic import is_list, is_pydantic_instance, is_curie_type, is_identifiable, is_literal
# from .models import class_refatt_map

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


def ga4gh_serialize2(vro: dict) -> str:
    """
    Provide a pure dict-based implementation of ga4gh_serialize using ga4gh_roll_up_refs
    """
    pass



def ga4gh_roll_up_refs(obj: dict, class_refatt_map={}):
    """
    Digestible key map is a map of unqualified class names to keys that can be rolled up.
    e.g. {"SequenceLocation": ["sequence"]}
    """
    obj_type = obj["type"]
    reffable_fields = class_refatt_map.get(obj_type, [])
    out = {}
    for k, v in obj.items():
        if k in reffable_fields:
            print(f"{k=}:{v=}")
            # three cases:
            # value is a ga4gh: CURIE, which can have digest extracted
            # value is a dict, which needs to itself be digested
            # value is something else
            if is_curie_type(v):
                pass



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
