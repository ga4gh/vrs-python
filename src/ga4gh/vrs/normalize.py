"""VRS object normalization functions

See https://vrs.ga4gh.org/en/stable/impl-guide/normalization.html

"""
import itertools
import logging
from enum import IntEnum
from typing import NamedTuple, Optional, Union

from bioutils.normalize import normalize as _normalize, NormalizationMode
from ga4gh.core import is_pydantic_instance, ga4gh_digest, pydantic_copy

from . import models
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
        len_trimmed_ref = len_trimmed_alt = 0
        # TODO: Return RLE for ref agree Alleles
    else:
        trim_ref_seq = ref_seq[trim_ival[0]: trim_ival[1]]
        trim_alt_seq = trim_alleles[1]
        len_trimmed_ref = len(trim_ref_seq)
        len_trimmed_alt = len(trim_alt_seq)

    # Compare the two allele sequences
    if not len_trimmed_ref and not len_trimmed_alt:
        return input_allele

    new_allele = pydantic_copy(input_allele)

    if len_trimmed_ref and len_trimmed_alt:
        new_allele.location.start = _get_new_allele_location_pos(
            trim_ival[0], start.pos_type
        )
        new_allele.location.end = _get_new_allele_location_pos(
            trim_ival[1], end.pos_type
        )
        new_allele.state.sequence = models.SequenceString(trim_alleles[1])
        return new_allele
    elif len_trimmed_ref:
        seed_length = len_trimmed_ref
    else:
        seed_length = len_trimmed_alt

    # Determine bounds of ambiguity
    new_ival, new_alleles = _normalize(
        ref_seq,
        trim_ival,
        (None, trim_alleles[1]),
        mode=NormalizationMode.EXPAND
    )

    new_allele.location.start = _get_new_allele_location_pos(
        new_ival[0], start.pos_type
    )
    new_allele.location.end = _get_new_allele_location_pos(
        new_ival[1], end.pos_type
    )

    extended_ref_seq = ref_seq[new_ival[0]: new_ival[1]]
    extended_alt_seq = new_alleles[1]

    if not extended_ref_seq:
        # If the reference sequence is empty this is an unambiguous insertion.
        # Return a new Allele with the trimmed alternate sequence as a Literal
        # Sequence Expression
        new_allele.state = models.LiteralSequenceExpression(
            sequence=models.SequenceString(extended_alt_seq)
        )
        return new_allele

    # Otherwise, determine if this is reference-derived (an RLE allele).
    len_extended_alt = len(extended_alt_seq)
    len_extended_ref = len(extended_ref_seq)

    if len_extended_alt < len_extended_ref:
        # If this is a deletion, it is reference-derived
        return _define_rle_allele(new_allele, len_extended_alt, seed_length, rle_seq_limit, extended_alt_seq)

    if len_extended_alt > len_extended_ref:
        # If this is an insertion, it may or may not be reference-derived.
        #
        # Determine the greatest factor `d` of the `seed length` such that `d`
        # is less than or equal to the length of the modified `reference sequence`,
        # and there exists a subsequence of length `d` derived from the modified
        # `reference sequence` that can be circularly expanded to recreate
        # the modified `alternate sequence`.
        factors = _factor_gen(seed_length)
        for cycle_length in factors:
            if cycle_length > len_extended_ref:
                continue
            cycle_start = len_extended_ref - cycle_length
            if _is_valid_cycle(cycle_start, extended_ref_seq, extended_alt_seq):
                return _define_rle_allele(
                    new_allele, len_extended_alt, cycle_length, rle_seq_limit, extended_alt_seq)

    new_allele.state = models.LiteralSequenceExpression(
        sequence=models.SequenceString(extended_alt_seq)
    )
    return new_allele


def _factor_gen(n):
    """Yields all factors of an integer `n`, in descending order"""
    lower_factors = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            yield n//i
            if n//i != i:
                lower_factors.append(i)
        i += 1
    for factor in reversed(lower_factors):
        yield factor


def _define_rle_allele(allele, length, repeat_subunit_length, rle_seq_limit, extended_alt_seq):
    # Otherwise, create the Allele as an RLE
    allele.state = models.ReferenceLengthExpression(
        length=length,
        repeatSubunitLength=repeat_subunit_length
    )

    if (rle_seq_limit and length <= rle_seq_limit) or (rle_seq_limit is None):
        allele.state.sequence = models.SequenceString(extended_alt_seq)

    return allele


def _is_valid_cycle(template_start, template, target):
    cycle = itertools.cycle(template[template_start:])
    for char in target[len(template):]:
        if char != next(cycle):
            return False
    return True

# TODO _normalize_genotype?


def _normalize_cis_phased_block(o, data_proxy=None):

    o.members = sorted(o.members, key=ga4gh_digest)
    return o


handlers = {
    "Allele": _normalize_allele,
    "CisPhasedBlock": _normalize_cis_phased_block,
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
