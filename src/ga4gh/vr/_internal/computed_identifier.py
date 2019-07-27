"""serializes, digests, and identifies GA4GH objects

"""

import logging
import re

import pkg_resources
import yaml

from ga4gh.core import sha512t24u
from .models import models

from canonicaljson import encode_canonical_json
import python_jsonschema_objects as pjs


__all__ = "ga4gh_digest ga4gh_identify ga4gh_serialize".split()

_logger = logging.getLogger(__name__)


# TODO: move to ga4gh.core.identfiers
schema_dir = pkg_resources.resource_filename(__name__, "data/schema")
cfg = yaml.safe_load(open(schema_dir + "/ga4gh.yaml"))
type_prefix_map = cfg["identifiers"]["type_prefix_map"]
namespace = cfg["identifiers"]["namespace"]
curie_sep = cfg["identifiers"]["curie_sep"]
ref_sep = cfg["identifiers"]["ref_sep"]
ga4gh_ir_regexp = re.compile(cfg["identifiers"]["regexp"])

ns_w_sep = namespace + curie_sep


def ga4gh_identify(vro):
    """return the GA4GH digest-based id for the object, as a CURIE
    (string)

    >>> import ga4gh.vr
    >>> interval = ga4gh.vr.models.Interval(start=10,end=11)
    >>> location = ga4gh.vr.models.Location(sequence_id="ga4gh:SQ.0123abcd", interval=interval)

    # Compute computed id: 
    >>> cid = ga4gh_identify(location)
    >>> cid
    ga4gh:SL.lGBsEujtdjKPTxBMCPAeGArLgNEuPN99

    """

    pfx = type_prefix_map[vro.type]
    digest = ga4gh_digest(vro)
    ir = f"{namespace}{curie_sep}{pfx}{ref_sep}{digest}"
    return ir


def ga4gh_digest(vro):
    """return the GA4GH digest for the object

    >>> import ga4gh.vr
    >>> interval = ga4gh.vr.models.Interval(start=10,end=11)
    >>> location = ga4gh.vr.models.Location(sequence_id="ga4gh:SQ.0123abcd", interval=interval)

    >>> ga4gh_digest(location)
    lGBsEujtdjKPTxBMCPAeGArLgNEuPN99

    """

    assert is_identifiable(vro), "ga4gh_digest called with non-identifiable object"
    if vro._digest is None:
        vro._digest = sha512t24u(ga4gh_serialize(vro))
    return vro._digest._value
    

def ga4gh_serialize(vro):
    """serialize object into a canonical format

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

    >>> import ga4gh.vr
    >>> interval = ga4gh.vr.models.Interval(start=10,end=11)
    >>> location = ga4gh.vr.models.Location(sequence_id="ga4gh:SQ.0123abcd", interval=interval)

    >>> ga4gh_serialize(location)
    b'{"interval":{"end":11,"start":10,"type":"SimpleInterval"},"sequence_id":"ga4gh:SQ.0123abcd","type":"SequenceLocation"}'

    """

    def dictify(vro, enref=True):
        """recursively converts (any) object to dictionary prior to
        serialization

        enref: if True, replace nested identifiable objects with
        digests ("enref" is opposite of "de-ref")

        """

        if vro is None:
            return None
        if is_literal(vro):
            v = vro._value
            if is_ga4gh_identifier(v):
                # strip ga4gh identifier to just the digest
                v = v.split(ref_sep, 1)[1]
            return v
        if is_class(vro):
            if is_identifiable(vro) and enref:
                return ga4gh_digest(vro)
            return {k: dictify(vro[k], enref=True)
                    for k in vro
                    if not (k.startswith("_") or vro[k] is None)}
        return vro


    # The canonicaljson package does everything we want. Use that with
    # the hope that it will be upward compatible with a future
    # ratified proposal for json canonicalization.
    #
    # The following alternative does the same thing for our use case.
    # It's included here as an outline for anyone implementing in
    # another language.  (canonicaljson escapes unicode characters, as
    # required by the VR spec, but this doesn't apply to any known
    # uses so these are equivalent.)

    # >> import json
    # >> def cjdump(a):
    # >>     return json.dumps(a, sort_keys=True, separators=(',',':'),
    #                          indent=None).encode("utf-8")

    vro_dict = dictify(vro, enref=False)
    return encode_canonical_json(vro_dict)


############################################################################
## INTERNAL
# TODO: Move these to utils in ga4gh.vr or perhaps ga4gh.core

def is_literal(vro):
    return isinstance(vro, pjs.literals.LiteralValue)

def is_class(vro):
    return isinstance(vro, pjs.classbuilder.ProtocolBase)

def is_identifiable(vro):
    return is_class(vro) and ("_digest" in vro)

def is_ga4gh_identifier(ir):
    return str(ir).startswith(ns_w_sep)

def parse_ga4gh_identifier(ir):
    return ga4gh_ir_regexp.match(str(ir)).groupdict()




if __name__ == "__main__":
    import ga4gh.vr
    interval = ga4gh.vr.models.Interval(start=10,end=11)
    location = ga4gh.vr.models.Location(sequence_id="ga4gh:SQ.0123abcd", interval=interval)
    state = ga4gh.vr.models.SequenceState(sequence="C")
    allele = ga4gh.vr.models.Allele(location=location, state=state)
    ga4gh_identify(location)
