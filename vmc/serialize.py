from ._const import enc, sep


def serialize(o):
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
        return "<{t}{sep}{o.namespace}{sep}{o.accession}>".format(sep=sep, t=t, o=o)

    if t == "Interval":
        return "<{t}{sep}{o.start}{sep}{o.end}>".format(sep=sep, t=t, o=o)

    if t == "Location":
        ival = serialize(o.interval)
        return "<{t}{sep}{o.sequence_id}{sep}{ival}>".format(sep=sep, t=t, o=o, ival=ival)

    if t == "Allele":
        return "<{t}{sep}{o.location_id}{sep}{o.state}>".format(sep=sep, t=t, o=o)

    if t == "Haplotype":
        # sort as well-defined binary encoding to circumvent locale-dependent sorting differences
        ids = sorted(i._value.encode(enc) for i in o.allele_ids)
        ids_str = ";".join(i.decode(enc) for i in ids)
        l = o.location_id or ""
        return "<{t}{sep}{l}{sep}{o.completeness}{sep}[{ids_str}]>".format(sep=sep, t=t, o=o, l=l, ids_str=ids_str)

    if t == "Genotype":
        # sort as well-defined binary encoding to circumvent locale-dependent sorting differences
        ids = sorted(i._value.encode(enc) for i in o.haplotype_ids)
        ids_str = ";".join(i.decode(enc) for i in ids)
        return "<{t}{sep}{o.completeness}{sep}[{ids_str}]>".format(sep=sep, t=t, o=o, ids_str=ids_str)

    raise Exception("Cannot serialize; unknown VMC object type: " + t)
