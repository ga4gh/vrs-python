"""serializes GA4GH VR objects

serialization converts a VR object (vro) into a *binary*
representation, typically in order to generate a digest.

Serialization implemented here occurs in two phases:
  [vro] ----dictify----> [dict] ----json-encoding----> [canonicaljson]


"""


import base64
import hashlib
import logging

from .const import ENC, NAMESPACE
from .models import models

from canonicaljson import encode_canonical_json
import python_jsonschema_objects as pjs


_logger = logging.getLogger(__name__)

ga4gh_model_prefixes = {
    # SQ: Sequence does not have a model
    models.Allele: "VA",
    models.Text: "VT",
    models.SequenceLocation: "SL",

    # The following are for post-1.0:
    # models.GeneLocation: "GL",
    # models.CytobandLocation: "CL",
    # models.Genotype: "VG",
    # models.Haplotype: "VH",
    # models.VariationSet: "VS",
}


def computed_id(o):
    """return the GA4GH digest-based id for the object, as a CURIE
    (string)

    >>> import ga4gh.vr
    >>> interval = ga4gh.vr.models.Interval(start=10,end=11)
    >>> location = ga4gh.vr.models.Location(sequence_id="GA4GH:GS_bogus", interval=interval)

    # Compute computed id: 
    >>> cid = computed_id(location)
    >>> cid
    'GA4GH:GL_RDaX1nGMg7D4M_Y9tiBQ_zG32cNkgkXQ'

    """

    if o.id is not None:
        return str(o.id)
    return "{ir.namespace}:{ir.accession}".format(ir=computed_identifier(o))


def computed_identifier(o):
    """return the GA4GH digest-based identifier for the object, as an Identifier

    >>> import ga4gh.vr
    >>> interval = ga4gh.vr.models.Interval(start=10,end=11)
    >>> location = ga4gh.vr.models.Location(sequence_id="GA4GH:GS_bogus", interval=interval)

    # Compute computed identifier: 
    >>> cid = computed_identifier(location)
    >>> cid
    <Identifier accession=GL_RDaX1nGMg7D4M_Y9tiBQ_zG32cNkgkXQ namespace=GA4GH>

    """

    pfx = ga4gh_model_prefixes[type(o)]
    gd = ga4gh_digest(serialize(o))
    ir = models.Identifier(namespace=NAMESPACE, accession=pfx + gd)
    return ir


def dictify(o, enref=True):
    """converts (any) object to dictionary prior to serialization

    enref: if True, replace nested identifiable objects with
    identifiers ("enref" is opposite of "de-ref")

    """

    def _dictify(o, enref, level):
        if o is None:
            return None
        if isinstance(o, pjs.literals.LiteralValue):
            return o._value
        if isinstance(o, pjs.classbuilder.ProtocolBase):
            if "id" in o and enref and level>0:
                # object is identifiable and user asked to enref
                if getattr(o, "id") is None:
                    setattr(o, "id", computed_id(o))
                return str(getattr(o, "id"))
            return {k: _dictify(o[k], enref, level+1) for k in o if o[k] is not None}
        _logger.critical(f"Got a {o} object")
        return o

    return _dictify(o=o, enref=enref, level=0)


def ga4gh_digest(blob):
    """generate a GA4GH digest for the given binary object

    A GA4GH digest is a convention for constructing and formatting
    digests for use as object identifiers. Specifically::
    
        * generate a SHA512 digest on binary data
        * truncate at 24 bytes
        * encode using base64url encoding

    Examples:
    >>> ga4gh_digest(b'')
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    >>> ga4gh_digest(b"ACGT")
    'aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'

    """

    digest_size = 24            # bytes
    digest = hashlib.sha512(blob).digest()
    tdigest_b64us = base64.urlsafe_b64encode(digest[:digest_size])
    return tdigest_b64us.decode("ASCII")


def serialize(o, enref=True):
    """serialize object using canonical json"""
    return encode_canonical_json(dictify(o, enref=enref))
