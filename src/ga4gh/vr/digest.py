"""Provides format and digest functions for VMC objects

A VMC digest is computed using a "truncated digest" (see below) on a
well-prescribed serialization of an object.

"""

import base64
import hashlib


from six import binary_type, text_type

from . import models
from ._const import ENC, NAMESPACE
from .serialize import serialize


abc = models                    # still needed?

vmc_model_prefixes = {
    # GS: Sequence does not have a model
    models.Allele: "GA",
    models.Genotype: "GG",
    models.Haplotype: "GH",
    models.SequenceLocation: "GL",
    models.GeneLocation: "GL",
    models.CytobandLocation: "GL",
    models.Text: "GT",
    models.VariationSet: "GVS",
}


def computed_id(o):
    """return the VMC digest-based id for the object, as a CURIE
    (string)

    >>> import ga4gh.vr
    >>> interval = ga4gh.vr.models.Interval(start=10,end=11)
    >>> location = ga4gh.vr.models.Location(sequence_id="VMC:GS_bogus", interval=interval)

    # Compute computed id: 
    >>> cid = computed_id(location)
    >>> cid
    'VMC:GL_RDaX1nGMg7D4M_Y9tiBQ_zG32cNkgkXQ'

    """

    if o.id is not None and o.id.startswith(NAMESPACE + ":"):
        return o.id

    return "{i.namespace}:{i.accession}".format(i=computed_identifier(o))


def computed_identifier(o):
    """return the VMC digest-based identifier for the object, as an Identifier

    >>> import ga4gh.vr
    >>> interval = ga4gh.vr.models.Interval(start=10,end=11)
    >>> location = ga4gh.vr.models.Location(sequence_id="VMC:GS_bogus", interval=interval)

    # Compute computed identifier: 
    >>> cid = computed_identifier(location)
    >>> cid
    <Identifier accession=GL_RDaX1nGMg7D4M_Y9tiBQ_zG32cNkgkXQ namespace=VMC>

    """

    pfx = vmc_model_prefixes.get(type(o), "??")
    dig = digest(o)
    accession = "{pfx}_{dig}".format(pfx=pfx, dig=dig)
    ir = models.Identifier(namespace=NAMESPACE, accession=accession)
    return ir


def digest(o, digest_size=24):
    """For a string object o, return the VMC digest of o.  For a VMC
    object o, return the VMC digest of the serialization of o.

    The VMC digest is an ASCII string derived by applying URL-safe
    Base64 encoding to the truncated SHA512 digest of the
    ASCII-encoded string.  Schematically, that is:
    
      digest := decode( b64us( trunc( sha512( encode(s) ), len )))

    or, in Python:

      digest = urlsafe_b64encode(sha512(s.encode("ascii")).digest()[:len]).decode()

    Example:
    >>> import ga4gh.vr

    >>> digest("")
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    >>> digest("", digest_size=12)
    'z4PhNX7vuL3xVChQ'

    >>> digest("ACGT")
    'aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'

    >>> interval = ga4gh.vr.models.Interval(start=10,end=11)
    >>> location = ga4gh.vr.models.Location(sequence_id="VMC:GS_bogus", interval=interval)

    >>> digest(serialize(location))
    'RDaX1nGMg7D4M_Y9tiBQ_zG32cNkgkXQ'

    >>> digest(location)
    'RDaX1nGMg7D4M_Y9tiBQ_zG32cNkgkXQ'

    """
    if isinstance(o, binary_type):
        b = o
    if isinstance(o, text_type):
        b = o.encode(ENC)
    else:
        b = serialize(o)
        if isinstance(b, text_type):
            b = b.encode(ENC)
    return _truncated_digest(b, digest_size=digest_size).decode(ENC)


############################################################################
# Internals

def _truncated_digest(blob, digest_size=24):
    """For the binary object blob, return the URL-safe, Base64 encoded,
    24-byte truncated SHA512 digest ascii encoded.  digest_size must
    be a multiple of 3.

    Examples:
    >>> _truncated_digest(b'')
    b'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    >>> _truncated_digest(b'', digest_size=24)
    b'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    >>> _truncated_digest(b'', digest_size=18)
    b'z4PhNX7vuL3xVChQ1m2AB9Yg'

    >>> _truncated_digest(b'', digest_size=9)
    b'z4PhNX7vuL3x'

    >>> _truncated_digest(b"ACGT", digest_size=24)
    b'aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'

    """

    # b64 encoding results in 4/3 size expansion of data and padded if
    # not multiple of 3, which doesn't make sense for this use
    assert digest_size % 3 == 0, "digest size must be multiple of 3"

    digest = hashlib.sha512(blob).digest()
    tdigest_b64us = base64.urlsafe_b64encode(digest[:digest_size])
    return tdigest_b64us
