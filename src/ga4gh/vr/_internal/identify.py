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

from .digest import ga4gh_digest
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



def identify(o):
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

    if o.id is not None:
        return str(o.id)
    pfx = _ga4gh_model_prefixes[type(o)]
    digest = ga4gh_digest(serialize(o))
    ir = f"{NAMESPACE}{PFX_REF_SEP}{pfx}{REF_SEP}{digest}"
    setattr(o, "id", ir)
    return ir


def serialize(o):
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

    return encode_canonical_json(_dictify(o))


############################################################################
## INTERNAL

def _dictify(o):
    """recursively converts (any) object to dictionary prior to
    serialization

    """

    def dictify_inner(o, enref=True):
        """enref: if True, replace nested identifiable objects with
        identifiers ("enref" is opposite of "de-ref")
        """
        if o is None:
            return None
        if isinstance(o, pjs.literals.LiteralValue):
            return o._value
        if isinstance(o, pjs.classbuilder.ProtocolBase):
            if "id" in o and enref:
                if o.id is None:
                    identify(o)
                return str(getattr(o, "id"))
            return {k: dictify_inner(o[k])
                    for k in o
                    if k != "id" and o[k] is not None}
        return o

    return dictify_inner(o=o, enref=False)  # don't enref first
