from . import models
from ._const import enc, namespace, sep


def serialize(o):
    """convert VMC object to canonical VMC serialized representation

    Serialization is the core of the VMC digest algorithm: Every VMC
    object must have a serialization in order for it to be
    identifiable (i.e., to be digested).

    Example:
    >>> import vmc
    >>> ir = vmc.models.Identifier(namespace="NCBI", accession="NC_000019.10")
    >>> serialize(ir)
    '<Identifier:NCBI:NC_000019.10>'

    >>> serialize("bogus")
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

    # TODO: this should be string id, not structured identifier
    if t == "Identifier":
        return "<{t}{sep}{o.namespace}{sep}{o.accession}>".format(sep=sep, t=t, o=o)

    if t == "Interval":
        return "<{t}{sep}{o.start}{sep}{o.end}>".format(sep=sep, t=t, o=o)

    if t == "Location":
        ir = _id_to_ir(o.sequence_id)
        return "<{t}{sep}{ir}{sep}{ival}>".format(sep=sep, t=t, o=o, ir=serialize(ir), ival=serialize(o.interval))

    if t == "Allele":
        ir = _id_to_ir(o.location_id)
        return "<{t}{sep}{ir}{sep}{o.state}>".format(sep=sep, t=t, ir=serialize(ir), o=o)

    if t == "Haplotype":
        # sort as well-defined binary encoding to circumvent locale-dependent sorting differences
        ids = sorted(serialize(_id_to_ir(str(i))).encode(enc) for i in o.allele_ids)
        l = o.location_id or ""
        return "<{t}{sep}{l}{sep}{o.completeness}:[{irss}]>".format(sep=sep, t=t, l=l, o=o, irss=";".join(i.decode(enc) for i in ids))

    if t == "Genotype":
        # sort as well-defined binary encoding to circumvent locale-dependent sorting differences
        ids = sorted(serialize(_id_to_ir(str(i))).encode(enc) for i in o.haplotype_ids)
        return "<{t}{sep}{o.completeness}:[{irss}]>".format(sep=sep, t=t, o=o, irss=";".join(i.decode(enc) for i in ids))

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
    assert id.startswith(namespace + ":"), "Must be a VMC id for this demo"
    ns, acc = id.split(":")
    return models.Identifier(namespace=ns, accession=acc)
