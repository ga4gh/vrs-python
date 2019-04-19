"""serializes VMC/GA4GH objects

Three serialization methods are implemented. The standard should
eventually choose only one, but there are complex tradeoffs left to
understand.  Choices:

1) Custom serialization
2) Canonical JSON
3) JSON Canonicalization Scheme

Briefly, 1 is under our control, but harder to keep up-to-date and
consistent with models. 2 Is implemented in, and pre-packaged, for
many languages but loosely spec'd. 3 is very well specified and on an
IETF approval track, but implemented by only one person so far and not
available pre-packaged. 

Reece's take: 3 is the most sound and, if IETF ratified, the best long
term option, but isn't available now. 1 is too painful. Go with 2 for
now, subject to verification across languages. Also, it will be
critical to have comprehensive tests of this functionality.

"""


import base64

from ._const import ENC, SEP

from canonicaljson import encode_canonical_json

def serialize_cj(o):
    """serialize object using canonical json"""
    # BUG: arrays need to be sorted for canonicalization
    # BUG: replace inline with ids (or vv, requires cache)
    d = o.as_dict()
    d.pop("id",None)
    return encode_canonical_json(d)


def serialize_vmc(o):
    """convert VMC object to canonical VMC serialized representation

    Serialization is the core of the VMC digest algorithm: Every VMC
    object must have a serialization in order for it to be
    identifiable (i.e., to be digested).

    Rules:
    * simple types (strings, ints) are serialized as-is
    * objects are serialized as <type|data>, where data might be a
      recursive serializations
    * pipe (|) is used to separate instance variables


    Example:
    >>> import vmc
    >>> iv = vmc.models.Interval(start=10, end=11)
    >>> serialize(iv)
    '<Interval|10|11>'

    >>> ir = vmc.models.Identifier(namespace="bob", accession="smith")
    >>> serialize(ir)
    '<Identifier|bob|smith>'

    >>> serialize("bob")
    Traceback (most recent call last):
    ...
    Exception: Cannot serialize; unknown VMC object type: str

    """

    # isinstance() fails here because nested classes built with
    # python_jsonschema_objects are coerced into the `abc` namespace.
    # So, we'll use the class "basename".
    t = o.__class__.__name__

    # Identifier serialization for completeness.
    if t == "Identifier":
        return "<{t}{sep}{o.namespace}{sep}{o.accession}>".format(sep=SEP, t=t, o=o)

    if t == "SimpleRegion":
        return "<{t}{sep}{o.start}{sep}{o.end}>".format(sep=SEP, t=t, o=o)

    if t == "SequenceLocation":
        rgn = serialize(o.region)
        return "<{t}{sep}{o.sequence_id}{sep}{rgn}>".format(sep=SEP, t=t, o=o, rgn=rgn)

    if t == "Text":
        b64 = base64.b64encode(str(o.definition).encode(ENC)).decode("ASCII")
        return f"<{t}{SEP}{b64}>"

    if t == "Allele":
        return "<{t}{sep}{o.location_id}{sep}{o.state}>".format(sep=SEP, t=t, o=o)

    if t == "Haplotype":
        # sort as well-defined binary encoding to circumvent locale-dependent sorting differences
        ids = sorted(i._value.encode(ENC) for i in o.allele_ids)
        ids_str = ";".join(i.decode(ENC) for i in ids)
        l = o.location_id or ""
        return "<{t}{sep}{l}{sep}{o.completeness}{sep}[{ids_str}]>".format(sep=SEP, t=t, o=o, l=l, ids_str=ids_str)

    if t == "Genotype":
        # sort as well-defined binary encoding to circumvent locale-dependent sorting differences
        ids = sorted(i._value.encode(ENC) for i in o.haplotype_ids)
        ids_str = ";".join(i.decode(ENC) for i in ids)
        return "<{t}{sep}{o.completeness}{sep}[{ids_str}]>".format(sep=SEP, t=t, o=o, ids_str=ids_str)

    raise Exception("Cannot serialize; unknown VMC object type: " + t)


serialize = serialize_vmc
