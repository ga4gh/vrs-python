"""VRS object normalization functions

See https://vrs.ga4gh.org/en/stable/impl-guide/normalization.html

"""
import logging
from enum import IntEnum
from itertools import cycle
from typing import NamedTuple, Optional, Union

from bioutils.normalize import normalize as _normalize, NormalizationMode
from ga4gh.core import is_pydantic_instance, ga4gh_digest, pydantic_copy

from ._internal import models
from .dataproxy import SequenceProxy


_logger = logging.getLogger(__name__)


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
    """Get a representative position for Alleles with Location start or end defined by Range

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
        pos_type = PosType.RANGE_LT_OR_EQUAL if pos0_is_none else PosType.RANGE_GT_OR_EQUAL

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


def _normalize_allele(input_allele, data_proxy, rle_seq_limit=50):
    """Normalize Allele using "fully-justified" normalization adapted from NCBI's
    VOCA. Fully-justified normalization expands such ambiguous representation over the
    entire region of ambiguity, resulting in an unambiguous representation that may be
    readily compared with other alleles.

    This function assumes that IRIs are dereferenced, providing a `SequenceReference` as
    the `allele.location.sequenceReference`. If a `SequenceReference` is not provided, the allele
    will be returned as is with no normalization.

    :param input_allele: Input VRS Allele object
    :param data_proxy: SeqRepo dataproxy
    :param rle_seq_limit: If RLE is set as the new state, set the limit for the length
        of the `sequence`.
        To exclude `sequence` from the response, set to 0.
        For no limit, set to `None`.

    Does not attempt to normalize Alleles with definite ranges and will instead return the
        `input_allele`
    """

    if isinstance(input_allele.location.sequenceReference, models.SequenceReference):
        alias = f"ga4gh:{input_allele.location.sequenceReference.refgetAccession}"
    else:
        _logger.warning(
            "`input_allele.location.sequenceReference` expects a `SequenceReference`, "
            "returning `input_allele` with no normalization."
        )
        return input_allele

    # Get reference sequence and interval
    ref_seq = SequenceProxy(data_proxy, alias)
    start = _get_allele_location_pos(input_allele, use_start=True)
    if start is None:
        return input_allele

    end = _get_allele_location_pos(input_allele, use_start=False)
    if end is None:
        return input_allele

    ival = (start.value, end.value)
    if input_allele.state.sequence:
        alleles = (None, input_allele.state.sequence.root)
    else:
        alleles = (None, "")

    # Trim common flanking sequence from Allele sequences.
    try:
        trim_ival, trim_alleles = _normalize(ref_seq, ival, alleles, mode=None, trim=True)
    except ValueError:
        # Occurs for ref agree Alleles (when alt = ref)
        len_ref_seq = len_alt_seq = 0
    else:
        trim_ref_seq = ref_seq[trim_ival[0]: trim_ival[1]]
        trim_alt_seq = trim_alleles[1]
        len_ref_seq = len(trim_ref_seq)
        len_alt_seq = len(trim_alt_seq)

    # Compare the two allele sequences
    if not len_ref_seq and not len_alt_seq:
        return input_allele

    new_allele = pydantic_copy(input_allele)

    if len_ref_seq and len_alt_seq:
        new_allele.location.start = _get_new_allele_location_pos(
            trim_ival[0], start.pos_type
        )
        new_allele.location.end = _get_new_allele_location_pos(
            trim_ival[1], end.pos_type
        )
        new_allele.state.sequence = models.SequenceString(trim_alleles[1])
        return new_allele

    # Determine bounds of ambiguity
    try:
        new_ival, new_alleles = _normalize(
            ref_seq,
            trim_ival,
            (None, trim_alleles[1]),
            mode=NormalizationMode.EXPAND
        )
    except ValueError:
        # Occurs for ref agree Alleles (when alt = ref)
        pass
    else:
        new_allele.location.start = _get_new_allele_location_pos(
            new_ival[0], start.pos_type
        )
        new_allele.location.end = _get_new_allele_location_pos(
            new_ival[1], end.pos_type
        )

        new_ref_seq = ref_seq[new_ival[0]: new_ival[1]]
        new_alt_seq = new_alleles[1]

        if not new_ref_seq or (
            new_alt_seq
            and _derive_seq_from_rle(
                ref_seq,
                new_allele.location.start,
                len_ref_seq or len_alt_seq,
                len(new_alt_seq),
            )
            != new_alt_seq
        ):
            # If the reference sequence is empty this is an unambiguous insertion or
            # If the reference sequence is not empty, but the alt sequence cannot be derived
            #   from the reference sequence
            # Return a new Allele with the trimmed alternate sequence as a Literal
            # Sequence Expression
            new_allele.state = models.LiteralSequenceExpression(
                sequence=models.SequenceString(new_alt_seq)
            )
        else:
            # Otherwise, return a new Allele using a RLE
            #  Use the ref/alt length from the trimmed allele, not the expanded form
            new_allele.state = models.ReferenceLengthExpression(
                length=len(new_alt_seq), repeatSubunitLength=len_ref_seq or len_alt_seq
            )

            if (rle_seq_limit and len(new_alt_seq) <= rle_seq_limit) or (
                rle_seq_limit is None
            ):
                new_allele.state.sequence = models.SequenceString(new_alt_seq)

    return new_allele


def _derive_seq_from_rle(
    ref_seq: SequenceProxy, start: int, repeatSubunitLength: int, length: int
):
    end = start + repeatSubunitLength
    subseq = ref_seq[start:end]
    c = cycle(subseq)
    derivedseq = ""
    for i in range(length):
        derivedseq += next(c)
    return derivedseq


# TODO _normalize_genotype?


def _normalize_haplotype(o, data_proxy=None):

    o.members = sorted(o.members, key=ga4gh_digest)
    return o


handlers = {
    "Allele": _normalize_allele,
    "Haplotype": _normalize_haplotype,
}


def normalize(vo, data_proxy=None, **kwargs):
    """normalize given vrs object, regardless of type

    kwargs:
        rle_seq_limit: If RLE is set as the new state, set the limit for the length
            of the `sequence`. To exclude `state.sequence`, set to 0.
    """

    assert is_pydantic_instance(vo)
    vo_type = vo.type

    if vo_type in handlers:
        handler = handlers[vo_type]
        return handler(vo, data_proxy, **kwargs)

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
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
            },
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "A",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }
    a = models.Allele(**allele_dict)

    allele2 = normalize(a, dp)

    a.state.sequence.root = "C"
    allele3 = normalize(a, dp)

    a.location.end = 44908823
    a.state.sequence.root = ""
    allele4 = normalize(a, dp)
