"""VRS object normalization functions

See https://vrs.ga4gh.org/en/stable/impl-guide/normalization.html

"""
from enum import IntEnum
from typing import NamedTuple, Optional, Union

from bioutils.normalize import normalize as _normalize, NormalizationMode
from ga4gh.core import is_pydantic_instance, ga4gh_digest, pydantic_copy

from ._internal import models
from .dataproxy import SequenceProxy


class PosType(IntEnum):
    """Define the kind of position on a location"""

    INTEGER = 0
    RANGE_LT_OR_EQUAL = 1
    RANGE_GT_OR_EQUAL = 2


class LocationPos(NamedTuple):
    """Define Allele location pos value and type"""

    value: int
    pos_type: PosType


def _get_allele_location_pos(
    allele_vo: models.Allele, use_start: bool = True
) -> Optional[LocationPos]:
    """Get Allele location start or end value for interval

    :param allele_vo: VRS Allele object
    :param use_start: `True` if using `allele_vo.location.start`. `False` if using
        `allele_vo.location.end`
    :return: A `LocationPos` if using integer or indefinite range. Otherwise return
        `None`
    """
    if use_start:
        pos = allele_vo.location.start
    else:
        pos = allele_vo.location.end

    if isinstance(pos, int):
        val = pos
        pos_type = PosType.INTEGER
    else:
        pos0_is_none = pos.root[0] is None
        pos1_is_none = pos.root[1] is None

        if not pos0_is_none and not pos1_is_none:  # definite range
            return None

        val = pos.root[0] or pos.root[1]
        if pos0_is_none:
            pos_type  = PosType.RANGE_LT_OR_EQUAL
        else:
            pos_type = PosType.RANGE_GT_OR_EQUAL

    return LocationPos(value=val, pos_type=pos_type)


def _get_new_allele_location_pos(
    new_pos_val: int, pos_type: PosType
) -> Union[int, models.Range]:
    """Get updated location pos on normalized allele

    :param new_pos_val: New position after normalization
    :param pos_type: Original position type used in pre-normalized VRS Allele object
    :return: Updated position as integer or VRS Range object
    """
    if pos_type == PosType.INTEGER:
        val = new_pos_val
    elif pos_type == PosType.RANGE_LT_OR_EQUAL:
        val = models.Range([None, new_pos_val])
    else:
        val = models.Range([new_pos_val, None])
    return val


def _normalize_allele(input_allele, data_proxy):
    """Normalize Allele using "fully-justified" normalization adapted from NCBI's
    VOCA. Fully-justified normalization expands such ambiguous representation over the
    entire region of ambiguity, resulting in an unambiguous representation that may be
    readily compared with other alleles.

    Does not attempt to normalize Allele's with definite ranges. Will return the
        `input_allele`
    """
    allele = pydantic_copy(input_allele)

    # Temporarily convert SequenceReference to IRI because it makes the code simpler.
    # This will be changed back to SequenceReference at the end of the method
    sequence_reference = None
    if isinstance(allele.location.sequence, models.SequenceReference):
        sequence_reference = allele.location.sequence
        allele.location.sequence = models.IRI(sequence_reference.refgetAccession)

    # Get reference sequence and interval
    ref_seq = SequenceProxy(data_proxy, allele.location.sequence.root)
    start = _get_allele_location_pos(allele, use_start=True)
    if start is None:
        return input_allele

    end = _get_allele_location_pos(allele, use_start=False)
    if end is None:
        return input_allele

    ival = (start.value, end.value)

    # Get alleles (the sequences to be normalized) for _normalize
    if allele.state.sequence:
        alleles = (None, allele.state.sequence.root)
    else:
        alleles = (None, "")

    # If one of Reference Allele Sequence or Alternate Allele Sequence is empty,
    # store the length of the non-empty sequence: this is the Repeat Subunit Length
    len_ref_seq = len(ref_seq[ival[0]: ival[1]])
    len_alt_seq = len(alleles[1])
    if not len_ref_seq and len_alt_seq:
        # Insertion
        repeat_subunit_len = len_alt_seq
    elif len_ref_seq and not len_alt_seq:
        # Deletion
        repeat_subunit_len = len_ref_seq
    else:
        repeat_subunit_len = 0

    new_allele = pydantic_copy(allele)
    try:
        new_ival, new_alleles = _normalize(ref_seq,
                                           ival,
                                           alleles=alleles,
                                           mode=NormalizationMode.EXPAND,
                                           anchor_length=0)

        new_allele.location.start = _get_new_allele_location_pos(
            new_ival[0], start.pos_type
        )
        new_allele.location.end = _get_new_allele_location_pos(
            new_ival[1], end.pos_type
        )
        new_ref_seq = ref_seq[new_ival[0]: new_ival[1]]

        if not new_ref_seq:
            # If the reference sequence is empty this is an unambiguous insertion.
            # Return a new Allele with the trimmed alternate sequence as a Literal
            # Sequence Expression
            new_allele.state = models.LiteralSequenceExpression(
                sequence=models.SequenceString(new_alleles[1])
            )
        else:
            # Otherwise, return a new Allele using a reference length expression, using
            # a Location specified by the coordinates of the new ival, a length
            # specified by the length of the alternate allele, and a repeat subunit
            # length
            allele.state = models.ReferenceLengthExpression(
                length=len(new_alleles[1]),
                sequence=models.SequenceString(new_ref_seq),
                repeat_subunit_length=repeat_subunit_len
            )
    except ValueError:
        # Occurs for ref agree Alleles (when alt = ref)
        pass

    # Convert IRI back to SequenceReference
    if sequence_reference:
        new_allele.location.sequence = sequence_reference

    return new_allele


# TODO _normalize_genotype?


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

    assert is_pydantic_instance(vo)
    vo_type = vo.type

    if vo_type in handlers:
        handler = handlers[vo_type]
        return handler(vo, data_proxy)

    # No handler for vo_type; pass-through unchanged
    return vo


if __name__ == "__main__":    # pragma: no cover
    # Requires seqrepo REST interface is running on this URL (e.g., using docker image)
    from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy
    seqrepo_rest_service_url = "http://localhost:5000/seqrepo"
    dp = SeqRepoRESTDataProxy(base_url=seqrepo_rest_service_url)

    # >>> dp.get_sequence("refseq:NC_000019.10", 44908820, 44908830)
    # " G C G C C T G G C A "
    #  |820      |825      | 830
    #
    allele_dict = {
        "location": {
            "end": 44908822,
            "start": 44908821,
            "sequence": "refseq:NC_000019.10",
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

    allele.location.interval.end = 44908823
    allele.state.sequence = ""
    allele4 = normalize(allele, dp)
