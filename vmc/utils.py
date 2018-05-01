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
