"""Uses SeqRepo (https://github.com/biocommons/biocommons.seqrepo) to
translate assigned identifiers (e.g., NC_000019.10 assigned by NCBI)
to VMC sequence digests (e.g.,
VMC:GS_IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl).

Although the VMC proposal requires sequence digests of the form
'VMC:GS_...', implementers are free to write the code to generate
these ids themselves.  The VMC demo uses SeqRepo because it includes
VMC digests in a lookup table.

"""

import os

import biocommons.seqrepo

import vmc


SEQREPO_ROOT_DIR = os.environ.get("SEQREPO_ROOT_DIR", "/usr/local/share/seqrepo")
SEQREPO_INSTANCE_NAME = os.environ.get("SEQREPO_INSTANCE", "testing")
seqrepo_instance_path = os.path.join(SEQREPO_ROOT_DIR, SEQREPO_INSTANCE_NAME)

_sr = None


def _get_seqrepo():
    global _sr
    if _sr is None:
        _sr = biocommons.seqrepo.SeqRepo(seqrepo_instance_path)
    return _sr


def get_vmc_sequence_identifier(identifier):
    """return VMC sequence Identifier (string) for a given Identifier from another namespace

    >>> get_vmc_sequence_identifier("RefSeq:NC_000019.10")
    'VMC:GS_IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl'

    >>> get_vmc_sequence_identifier("RefSeq:bogus")
    Traceback (most recent call last):
    ...
    KeyError: 'refseq:bogus'

    # also accepts an Identifier
    >>> from vmc import models
    >>> ir = models.Identifier(namespace="RefSeq", accession="NC_000019.10")
    >>> get_vmc_sequence_identifier(ir)
    'VMC:GS_IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl'

    """

    _sr = _get_seqrepo()
    if isinstance(identifier, vmc.models.Identifier):
        identifier = "{i.namespace}:{i.accession}".format(i=identifier)
    return _sr.translate_identifier(identifier, target_namespaces=["VMC"])[0]


def get_reference_sequence(id, start, end):
    return _get_seqrepo()[id][start:end]
