"""Provides format and digest functions for VMC objects

A VMC digest is computed using a "truncated digest" (see below) on a
well-prescribed serialization of an object.

"""

import base64
import hashlib

from . import models
from ._const import enc, namespace
from .serialize import serialize


abc = models                    # still needed?

vmc_model_prefixes = {
    # GS: Sequence does not have a model
    models.Allele: "GA",
    models.Genotype: "GG",
    models.Haplotype: "GH",
    models.Location: "GL",
}


def computed_id(o):
    """return the VMC computed identifier for the object, as a CURIE
    (string)

    """

    if o.id is not None and o.id.startswith(namespace + ":"):
        return o.id

    return "{i.namespace}:{i.accession}".format(i=computed_identifier(o))


def computed_identifier(o):
    """return the VMC computed identifier for the object, as an Identifier

    >>> import vmc
    >>> interval = vmc.models.Interval(start=10,end=11)
    >>> location = vmc.models.Location(sequence_id="VMC:GS_bogus", interval=interval)

    # Compute computed identifier: 
    >>> cid = computed_id(location)
    >>> cid
    'VMC:GL_LBe0tQtFb1wtDpiwDMM1ixBIYYn171fT'

    # Setting the id will preempt computing the identifier again:
    >>> location.id = cid
    >>> cid = computed_id(location)
    >>> cid
    'VMC:GL_LBe0tQtFb1wtDpiwDMM1ixBIYYn171fT'

    """

    pfx = vmc_model_prefixes[type(o)]
    dig = digest(o)
    accession = "{pfx}_{dig}".format(pfx=pfx, dig=dig)
    ir = models.Identifier(namespace=namespace, accession=accession)
    return ir


def digest(o):
    """For a VMC object o, return the URL-safe, Base64 encoded, 24-byte
    truncated SHA512 digest as unicode

    Example:
    >>> import vmc
    >>> ir = vmc.models.Identifier(namespace="NCBI", accession="NC_000019.10")
    >>> digest(ir)
    'Kh-Ml83IE-vt7Fu9XLBDmWWJlRzkOvcF'

    """
    ser = serialize(o)
    return _truncated_digest(ser.encode(enc)).decode(enc)


############################################################################
# Internals

def _truncated_digest(blob):
    """For the binary object blob, return the URL-safe, Base64 encoded,
    24-byte truncated SHA512 digest as unicode

    Example:
    >>> _truncated_digest(b'')
    b'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    """
    digest = hashlib.sha512(blob).digest()
    tdigest_b64us = base64.urlsafe_b64encode(digest[:24])
    return tdigest_b64us
