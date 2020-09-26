import logging

from bioutils.normalize import normalize as _normalize, NormalizationMode

from ..core import is_pjs_instance, pjs_copy
from ._internal import models
from .dataproxy import SequenceProxy


_logger = logging.getLogger(__name__)



def _normalize_allele(allele, data_proxy):
    sequence = SequenceProxy(data_proxy, allele.location.sequence_id._value)
    ival = (allele.location.interval.start._value, allele.location.interval.end._value)
    alleles = (None, allele.state.sequence._value)

    new_allele = pjs_copy(allele)

    try:
        new_ival, new_alleles = _normalize(sequence, ival,
                                           alleles=alleles,
                                           mode=NormalizationMode.EXPAND,
                                           anchor_length=0)
        new_allele.location.interval.start = new_ival[0]
        new_allele.location.interval.end = new_ival[1]
        new_allele.state.sequence = new_alleles[1]
    except ValueError:
        # Occurs for ref agree Alleles (when alt = ref)
        pass

    return new_allele


handlers = {
    "Allele": _normalize_allele,
}


def normalize(vo, data_proxy):
    assert is_pjs_instance(vo)

    vo_type = vo.type._value

    if vo_type in handlers:
        handler = handlers[vo_type]
        return handler(vo, data_proxy)
    else:
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
        'location': {
            'interval': {
                'end': 44908822,
                'start': 44908821,
                'type': 'SimpleInterval'
            },
            'sequence_id': 'refseq:NC_000019.10',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'A',
            'type': 'SequenceState'
        },
        'type': 'Allele'
    }
    allele = models.Allele(**allele_dict)
    

    allele2 = normalize(allele, dp)
    
    allele.state.sequence = "C"
    allele3 = normalize(allele, dp)
    
    allele.location.interval.end = 44908823
    allele.state.sequence = ""
    allele4 = normalize(allele, dp)
