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
