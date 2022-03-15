"""VRS object normalization functions

See https://vrs.ga4gh.org/en/stable/impl-guide/normalization.html

"""


import logging

from bioutils.normalize import normalize as _normalize, NormalizationMode
from ga4gh.core import is_pjs_instance, pjs_copy, ga4gh_digest
from python_jsonschema_objects.validators import ValidationError

from ._internal import models
from .dataproxy import SequenceProxy


_logger = logging.getLogger(__name__)


def _normalize_allele(allele, data_proxy):
    sequence = SequenceProxy(data_proxy, allele.location.sequence_id._value)

    _interval_type = allele.location.interval.type
    if _interval_type == "SimpleInterval":
        ival = (allele.location.interval.start._value, allele.location.interval.end._value)
    elif _interval_type == "SequenceInterval":
        ival = (allele.location.interval.start.value, allele.location.interval.end.value)

    _allele_state = allele.state.type
    _states_with_sequence = ['SequenceState', 'LiteralSequenceExpression']
    if _allele_state in _states_with_sequence:
        alleles = (None, allele.state.sequence._value)
    elif _allele_state == 'RepeatedSequenceExpression' and \
            allele.state.seq_expr.type in _states_with_sequence:
        alleles = (None, allele.state.seq_expr.sequence._value)
    else:
        alleles = (None, '')

    new_allele = pjs_copy(allele)

    try:
        new_ival, new_alleles = _normalize(sequence, ival,
                                           alleles=alleles,
                                           mode=NormalizationMode.EXPAND,
                                           anchor_length=0)

        _new_allele_ival_type = new_allele.location.interval.type
        if _new_allele_ival_type == "SimpleInterval":
            new_allele.location.interval.start = new_ival[0]
            new_allele.location.interval.end = new_ival[1]
        elif _new_allele_ival_type == "SequenceInterval":
            new_allele.location.interval.start.value = new_ival[0]
            new_allele.location.interval.end.value = new_ival[1]

        if new_allele.state.type in _states_with_sequence:
            new_allele.state.sequence = new_alleles[1]
    except ValueError:
        # Occurs for ref agree Alleles (when alt = ref)
        pass

    return new_allele

def _normalize_haplotype(o, data_proxy=None):
    o.members = sorted(o.members, key=ga4gh_digest)
    return o

def _normalize_variationset(o, data_proxy=None):
    o.members = sorted(o.members, key=ga4gh_digest)
    return o


handlers = {
    "Allele": _normalize_allele,
    "Haplotype": _normalize_haplotype,
    "VariationSet": _normalize_variationset,
}


def normalize(vo, data_proxy=None):
    """normalize given vrs object, regardless of type"""

    assert is_pjs_instance(vo)
    vo_type = vo.type._value

    if vo_type in handlers:
        handler = handlers[vo_type]
        return handler(vo, data_proxy)

    # No handler for vo_type; pass-through unchanged
    return vo



if __name__ == "__main__":      # pragma: no cover
    # Requires seqrepo REST interface is running on this URL (e.g., using docker image)
    from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy
    seqrepo_rest_service_url = "http://localhost:5000/seqrepo"
    dp = SeqRepoRESTDataProxy(base_url=seqrepo_rest_service_url)


    # >>> dp.get_sequence("refseq:NC_000019.10", 44908820, 44908830)
    # ' G C G C C T G G C A '
    #  |820      |825      | 830
    #
    allele_dict = {
        "location": {
            "interval": {
                "end": {'value': 44908822, 'type': 'Number'},
                "start": {'value': 44908821, 'type': 'Number'},
                "type": "SequenceInterval"
            },
            "sequence_id": "refseq:NC_000019.10",
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "A",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }
    allele = models.Allele(**allele_dict)


    allele2 = normalize(allele, dp)

    allele.state.sequence = "C"
    allele3 = normalize(allele, dp)

    allele.location.interval.end.value = 44908823
    allele.state.sequence = ""
    allele4 = normalize(allele, dp)
