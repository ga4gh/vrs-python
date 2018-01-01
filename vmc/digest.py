"""Provides format and digest functions for VMC objects

A VMC digest is computed using a "truncated digest" (see below) on a
well-prescribed serialization of an object.

"""

import base64
import hashlib

from . import models

enc = "ASCII"

abc = models                    # still needed?

vmc_namespace = "VMC"
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

    if o.id is not None and o.id.startswith(vmc_namespace + ":"):
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
    dig = vmc_digest(o)
    accession = "{pfx}_{dig}".format(pfx=pfx, dig=dig)
    ir = models.Identifier(namespace=vmc_namespace, accession=accession)
    return ir


def vmc_digest(o):
    """For a VMC object o, return the URL-safe, Base64 encoded, 24-byte
    truncated SHA512 digest as unicode

    Example:
    >>> import vmc
    >>> ir = vmc.models.Identifier(namespace="NCBI", accession="NC_000019.10")
    >>> vmc_digest(ir)
    'Kh-Ml83IE-vt7Fu9XLBDmWWJlRzkOvcF'

    """
    ser = vmc_serialize(o)
    return _truncated_digest(ser.encode(enc)).decode(enc)


def vmc_serialize(o):
    """convert VMC object to canonical VMC serialized representation

    Serialization is the core of the VMC digest algorithm: Every VMC
    object must have a serialization in order for it to be
    identifiable (i.e., to be digested).

    Example:
    >>> import vmc
    >>> ir = vmc.models.Identifier(namespace="NCBI", accession="NC_000019.10")
    >>> vmc_serialize(ir)
    '<Identifier:NCBI:NC_000019.10>'

    >>> vmc_serialize("bogus")
    Traceback (most recent call last):
    ...
    Exception: Unknown type: str

    Two wrinkles:

    * The digest must be based on VMC (external) identifiers. For this
      demo, the VMC identifiers and internal ids are equivalent, but
      that need not be true generally.  For example, an implementation
      might choose to use uuids for ids, with VMC identifiers
      attached.  (See _id_to_ir for shortcut used in this demo.)

    * Alleles need to be reliably ordered.  Unfortunately, URL-safe
      Base 64 encodings use characters that are subject to
      locale-dependent sort differences.  So, we sort here by the
      ASCII-encoded form (essentially, binary sorted).

    """

    # isinstance() fails here because nested classes built with
    # python_jsonschema_objects are coerced into the `abc` namespace.
    # So, we'll use the class "basename".
    t = o.__class__.__name__

    if t == "Identifier":
        return "<{t}:{o.namespace}:{o.accession}>".format(t=t, o=o)

    if t == "Interval":
        return "<{t}:{o.start}:{o.end}>".format(t=t, o=o)

    if t == "Location":
        ir = _id_to_ir(o.sequence_id)
        return "<{t}:{ir}:{ival}>".format(t=t, ir=vmc_serialize(ir), ival=vmc_serialize(o.interval))

    if t == "Allele":
        ir = _id_to_ir(o.location_id)
        return "<{t}:{ir}:{o.state}>".format(t=t, ir=vmc_serialize(ir), o=o)

    if t == "Haplotype":
        # sort as well-defined binary encoding to circumvent locale-dependent sorting differences
        ids = sorted(vmc_serialize(_id_to_ir(str(i))).encode(enc) for i in o.allele_ids)
        return "<{t}:{o.completeness}:[{irss}]>".format(t=t, o=o, irss=";".join(i.decode(enc) for i in ids))

    if t == "Genotype":
        # sort as well-defined binary encoding to circumvent locale-dependent sorting differences
        ids = sorted(vmc_serialize(_id_to_ir(str(i))).encode(enc) for i in o.haplotype_ids)
        return "<{t}:{o.completeness}:[{irss}]>".format(t=t, o=o, irss=";".join(i.decode(enc) for i in ids))

    raise Exception("Unknown type: " + t)


############################################################################
# Internals

def _id_to_ir(id):
    """Convert internal id to identifier string.

    For the VMC demo, we're going to assume that the computed
    identifier is used as the id, and therefore we can just assert
    that the id begins with the "VMC:" and then return a synthesized
    Identifier.

    If an implementation uses a different internal id, this function
    needs to be changed.

    """
    assert id.startswith(vmc_namespace + ":"), "Must be a VMC id for this demo"
    ns, acc = id.split(":")
    return models.Identifier(namespace=ns, accession=acc)


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
