"""serializes GA4GH VR objects

Serialization implemented here occurs in these phases:

   [vro] -----> [dict] -----> [cjson] -----> [digest] -----> [computed_id]
     |--dictify---|--encode_cj---|-ga4gh_digest-|------<format id>-------|
     |--------serialize----------|
     |----------------------------identify-------------------------------|

[vro] = VR object
[dict] = Python dict
[cjson] = canonical JSON format
[digest] = ga4gh_digest result
[computed_id] = CURIE-formatted string identifier

dictify: convert vro to dict, by default replacing inlined objects
with identifiers ("enref" option)

encode_cj: encode dict as UTF-8 encoded JSON per spec

serialization: dictify + encode_cj; converts a VR object (vro) into a
*binary* representation, typically in order to generate a digest.

"""


import logging

from ga4gh.core import ga4gh_digest
from .models import models

from canonicaljson import encode_canonical_json
import python_jsonschema_objects as pjs


__all__ = "identify ga4gh_digest serialize".split()


_logger = logging.getLogger(__name__)

_ga4gh_model_prefixes = {
    # SQ: Sequence does not have a model
    models.Allele: "VA",
    models.SequenceLocation: "SL",
    models.Text: "VT",

    # The following are for post-1.0:
    # models.CytobandLocation: "CL",
    # models.GeneLocation: "GL",
    # models.Genotype: "VG",
    # models.Haplotype: "VH",
    # models.VariationSet: "VS",
}

# computed identifer format:
# <NAMESPACE><PFX_REF_SEP><type-prefix><REF_SEP><digest>
# eg., ga4gh:SQ/0123abcd
NAMESPACE = "ga4gh"
PFX_REF_SEP = ":"
REF_SEP = "/"


def identify(vro):
    """return the GA4GH digest-based id for the object, as a CURIE
    (string)

    >>> import ga4gh.vr
    >>> interval = ga4gh.vr.models.Interval(start=10,end=11)
    >>> location = ga4gh.vr.models.Location(sequence_id="GA4GH:GS_bogus", interval=interval)

    # Compute computed id: 
    >>> cid = identify(location)
    >>> cid
    'GA4GH:GLRDaX1nGMg7D4M_Y9tiBQ_zG32cNkgkXQ'

    """

    if vro.id is not None:
        return str(vro.id)
    pfx = _ga4gh_model_prefixes[type(vro)]
    digest = ga4gh_digest(serialize(vro))
    ir = f"{NAMESPACE}{PFX_REF_SEP}{pfx}{REF_SEP}{digest}"
    setattr(vro, "id", ir)
    return ir


def serialize(vro):
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

    See https://github.com/ga4gh/vr-schema/issues/31 for a discussion
    """

    # The canonicaljson pacakge does everything we want. Use that with
    # the expectation (hope) that it will be upward compatible with a
    # future ratified proposal.
    #
    # The following alternative seems to do the same thing for our use
    # case.  It's included here as an outline for anyone implementing
    # in another language.
    # >> import json
    # >> def cjdump(a):
    # >>     return json.dumps(a, sort_keys=True, separators=(',',':'),
    #                          indent=None).encode("utf-8")

    return encode_canonical_json(_dictify(vro))


############################################################################
## INTERNAL


def is_literal(vro):
    return isinstance(vro, pjs.literals.LiteralValue)

def is_class(vro):
    return isinstance(vro, pjs.classbuilder.ProtocolBase)

def is_identifiable(vro):
    return is_class(vro) and ("id" in vro)

def _dictify(vro):
    """recursively converts (any) object to dictionary prior to
    serialization

    """

    def dictify_inner(vro, enref=True):
        """enref: if True, replace nested identifiable objects with
        identifiers ("enref" is opposite of "de-ref")
        """
        if vro is None:
            return None
        if is_literal(vro):
            return vro._value
        if is_class(vro):
            if "id" in vro and enref:
                if vro.id is None:
                    identify(vro)
                elif not str(vro.id).startswith(NAMESPACE + PFX_REF_SEP):
                    raise GA4GHError("ga4gh computed identifiers require that nested objects use GA4GH computed identifiers") 
                return str(getattr(vro, "id"))
            return {k: dictify_inner(vro[k])
                    for k in vro
                    if k != "id" and vro[k] is not None}
        return vro

    return dictify_inner(vro=vro, enref=False)  # don't enref first
