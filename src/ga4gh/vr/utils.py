import json

import python_jsonschema_objects
from bioutils.accessions import infer_namespace

from ._models import models



def id_to_ir(id):
    """Convert internal id to identifier string.

    For the VMC demo, we're going to assume that the computed
    identifier is used as the id, and therefore we can just assert
    that the id begins with the "VMC:" and then return a synthesized
    Identifier.

    If an implementation uses a different internal id, this function
    needs to be changed.

    >>> id_to_ir("VMC:bogus")
    <Identifier accession=bogus namespace=VMC>

    """
    ns, acc = id.split(":")
    return models.Identifier(namespace=ns, accession=acc)


def ir_to_id(ir):
    return "{ir.namespace}:{ir.accession}".format(ir=ir)


def json_pretty_format(j):
    """pretty print object as json"""
    return json.dumps(json.loads(j), indent=4, sort_keys=True, ensure_ascii=False)


def is_vr_instance(o):
    return isinstance(o, python_jsonschema_objects.classbuilder.ProtocolBase)


def coerce_namespace(ac):
    """given an accession, prefix with inferred namespace if not present

    >>> coerce_namespace("refseq:NM_01234.5")
    'refseq:NM_01234.5'

    >>> coerce_namespace("NM_01234.5")
    'refseq:NM_01234.5'

    >>> coerce_namespace("QQ_01234.5")
    Traceback (most recent call last):
    ...
    ValueError: Could not infer namespace for QQ_01234.5

    >>> coerce_namespace("bogus:QQ_01234.5")
    'bogus:QQ_01234.5'

    """
    if ":" not in ac:
        ns = infer_namespace(ac)
        if ns is None:
            raise ValueError(f"Could not infer namespace for {ac}")
        ac = ns + ":" + ac
    return ac
