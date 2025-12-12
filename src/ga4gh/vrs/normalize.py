"""VRS object normalization functions

See https://vrs.ga4gh.org/en/stable/impl-guide/normalization.html

"""

import itertools
import logging
from enum import IntEnum
from typing import NamedTuple

from bioutils.normalize import NormalizationMode
from bioutils.normalize import normalize as _normalize

from ga4gh.core import ga4gh_digest, is_pydantic_instance, pydantic_copy
from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import SequenceProxy, _DataProxy

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
) -> LocationPos | None:
    """Get a representative position for Alleles with Location start or end defined by Range

    :param allele_vo: VRS Allele object
    :param use_start: `True` if using `allele_vo.location.start`. `False` if using
        `allele_vo.location.end`
    :return: A `LocationPos` if using integer or indefinite range. Otherwise return
        `None`
    """
    pos = allele_vo.location.start if use_start else allele_vo.location.end

    if isinstance(pos, int):
        val = pos
        pos_type = PosType.INTEGER
    else:
        pos0_is_none = pos.root[0] is None
        pos1_is_none = pos.root[1] is None

        if not pos0_is_none and not pos1_is_none:  # definite range
            return None

        val = pos.root[0] or pos.root[1]
        pos_type = (
            PosType.RANGE_LT_OR_EQUAL if pos0_is_none else PosType.RANGE_GT_OR_EQUAL
        )

    return LocationPos(value=val, pos_type=pos_type)


def _get_new_allele_location_pos(
    new_pos_val: int, pos_type: PosType
) -> int | models.Range:
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


def _normalize_allele(input_allele: models.Allele, data_proxy, rle_seq_limit=50):
    """Normalize Allele using "fully-justified" normalization adapted from NCBI's
    VOCA. Fully-justified normalization expands such ambiguous representation over the
    entire region of ambiguity, resulting in an unambiguous representation that may be
    readily compared with other alleles.

    This function assumes that IRIs are dereferenced, providing a `SequenceReference` as
    the `allele.location.sequenceReference`. If a `SequenceReference` is not provided, the allele
    will be returned as is with no normalization.

    Does not attempt to normalize Alleles with definite ranges or non-LiteralSequenceExpression
    states and will instead return the `input_allele`.

    If attempting to re-normalize a ReferenceLengthExpression allele,
    use denormalize_reference_length_expression to construct a LiteralSequenceExpression
    allele and then call normalize() on that LiteralSequenceExpression allele.

    See: https://vrs.ga4gh.org/en/2.0/conventions/normalization.html#literalsequenceexpression-alleles

    :param input_allele: Input VRS Allele object
    :param data_proxy: SeqRepo dataproxy
    :param rle_seq_limit: If RLE is set as the new state, set the limit for the length
        of the `sequence`.
        To exclude `sequence` from the response, set to 0.
        For no limit, set to `None`.
    """
    # Algorithm applies to LiteralSequenceExpression alleles only; other states are returned unchanged
    if not isinstance(input_allele.state, models.LiteralSequenceExpression):
        _logger.warning(
            "`input_allele.state` was not a LiteralSequenceExpression, returning `input_allele` with no normalization."
        )
        return input_allele

    if isinstance(input_allele.location.sequenceReference, models.SequenceReference):
        alias = f"ga4gh:{input_allele.location.sequenceReference.refgetAccession}"
    else:
        _logger.warning(
            "`input_allele.location.sequenceReference` expects a `SequenceReference`, "
            "returning `input_allele` with no normalization."
        )
        return input_allele

    # 0: Get reference sequence and interval
    ref_seq = SequenceProxy(data_proxy, alias)
    start = _get_allele_location_pos(input_allele, use_start=True)
    if start is None:
        return input_allele

    end = _get_allele_location_pos(input_allele, use_start=False)
    if end is None:
        return input_allele

    ival = (start.value, end.value)
    start_pos_type = start.pos_type
    end_pos_type = end.pos_type
    alt_seq = input_allele.state.sequence.root or ""
    alleles = (None, alt_seq)

    # 1: trim shared flanking sequence
    trim_ival, trim_alleles = _trim_for_normalization(
        ref_seq, ival, alleles, start, end
    )

    trim_ref_seq = ref_seq[trim_ival[0] : trim_ival[1]]
    trim_alt_seq = trim_alleles[1]
    len_trimmed_ref = len(trim_ref_seq)
    len_trimmed_alt = len(trim_alt_seq)
    seed_length = len_trimmed_ref or len_trimmed_alt
    identity_case = trim_ref_seq == trim_alt_seq

    new_allele: models.Allele = pydantic_copy(input_allele)

    # 2.a: Reference allele (ref==alt after trim): use original span and return RLE
    # length = repeatSubunitLength = seed_length (the input sequence length)
    if identity_case:
        _set_location_from_interval(new_allele, ival, start_pos_type, end_pos_type)
        return _define_rle_allele(
            new_allele, seed_length, seed_length, rle_seq_limit, alt_seq
        )

    # 2.b: Substitution: both sides non-empty and different after trim.
    if len_trimmed_ref and len_trimmed_alt:
        _set_location_from_interval(new_allele, trim_ival, start_pos_type, end_pos_type)
        new_allele.state.sequence = models.sequenceString(trim_alt_seq)
        return new_allele

    # 3: Expand ambiguity by rolling left + right
    new_ival, new_alleles = _normalize(
        ref_seq,
        trim_ival,
        (None, trim_alt_seq),
        mode=NormalizationMode.EXPAND,
        trim=not identity_case,  # bioutils will not trim identical alleles
    )

    # 4: Get the extended sequences
    extended_ref_seq = ref_seq[new_ival[0] : new_ival[1]]
    extended_alt_seq = new_alleles[1]
    len_extended_alt = len(extended_alt_seq)
    len_extended_ref = len(extended_ref_seq)

    # 5.a
    if not extended_ref_seq:
        _set_location_from_interval(new_allele, new_ival, start_pos_type, end_pos_type)
        new_allele.state = models.LiteralSequenceExpression(
            sequence=models.sequenceString(extended_alt_seq)
        )
        return new_allele

    # 5.b
    if len_extended_alt < len_extended_ref:
        _set_location_from_interval(new_allele, new_ival, start_pos_type, end_pos_type)
        return _define_rle_allele(
            new_allele, len_extended_alt, seed_length, rle_seq_limit, extended_alt_seq
        )

    # 5.c
    if len_extended_alt > len_extended_ref:
        factors = _factor_gen(seed_length)
        for cycle_length in factors:
            if cycle_length > len_extended_ref:
                continue
            cycle_start = len_extended_ref - cycle_length
            if _is_valid_cycle(cycle_start, extended_ref_seq, extended_alt_seq):
                _set_location_from_interval(
                    new_allele, new_ival, start_pos_type, end_pos_type
                )
                # 5.c.2 / 5.d: reference-derived ambiguous insertion
                return _define_rle_allele(
                    new_allele,
                    len_extended_alt,
                    cycle_length,
                    rle_seq_limit,
                    extended_alt_seq,
                )
        # 5.c.3
        _set_location_from_interval(new_allele, new_ival, start_pos_type, end_pos_type)
        new_allele.state = models.LiteralSequenceExpression(
            sequence=models.sequenceString(extended_alt_seq)
        )
        return new_allele

    # 5.e: Otherwise return literal Allele using expanded interval/state (spec step 5 final bullet)
    _set_location_from_interval(new_allele, new_ival, start_pos_type, end_pos_type)
    new_allele.state = models.LiteralSequenceExpression(
        sequence=models.sequenceString(extended_alt_seq)
    )
    return new_allele


def _trim_for_normalization(ref_seq: SequenceProxy, ival: tuple[int, int], alleles: tuple[None, str], start: LocationPos, end: LocationPos):
    """Trim common prefix and suffix from the intervals.

    Return the trimmed interval and trimmed alleles:
        ((trim_start, trim_end), (trim_ref, trim_alt))

    If the alleles are the same, return the original interval and alleles.
    The first allele (ref) will be populated by bioutils.
    """
    try:
        trim_ival, trim_alleles = _normalize(
            ref_seq, ival, alleles, mode=None, trim=True
        )
    except ValueError as e:
        ref_at_location = ref_seq[start.value : end.value]
        alt_seq = alleles[1]
        if ref_at_location == alt_seq:
            # return (ival, (None, ref_at_location))
            return ival, alleles
        msg = (
            "Unexpected bioutils trim error for non reference allele: "
            f"ref='{ref_at_location}', alt='{alt_seq}'"
        )
        raise ValueError(msg) from e

    return trim_ival, trim_alleles


def _set_location_from_interval(allele: models.Allele, ival: tuple[int, int], start_pos_type: PosType, end_pos_type: PosType) -> None:
"""Update ``allele`` start and end location""" 
    allele.location.start = _get_new_allele_location_pos(ival[0], start_pos_type)
    allele.location.end = _get_new_allele_location_pos(ival[1], end_pos_type)


def denormalize_reference_length_expression(
    ref_seq: str,
    repeat_subunit_length: int,
    alt_length: int,
) -> str:
    """Reverse the process of repeat subunit compaction for a ReferenceLengthExpression.
    The alt is reference-derived so repeat the repeat subunit until it is it to alt_length
    characters long. May include a trailing partial copy of the repeat subunit.

    e.g. "ACGT" (repeat_subunit_length=4, length=10) -> "ACGTACGTAC" (ACGT ACGT AC)

    :param ref_seq: The reference sequence to be denormalized.
    :param repeat_subunit_length: The length of the repeat subunit.
    :param alt_length: The length of the alternate sequence that was compacted during normalization.
    """
    repeat_subunit = ref_seq[:repeat_subunit_length]
    if len(repeat_subunit) != repeat_subunit_length:
        raise ValueError(  # noqa: TRY003
            f"Repeat subunit length {repeat_subunit_length} is not equal to the length of the repeat subunit {len(repeat_subunit)}"  # noqa: EM102
        )
    repeat_count = alt_length // repeat_subunit_length
    remainder = alt_length % repeat_subunit_length
    alt = repeat_subunit * repeat_count
    alt += repeat_subunit[:remainder]
    return alt


def _factor_gen(n):
    """Yield all factors of an integer `n`, in descending order"""
    lower_factors = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            yield n // i
            if n // i != i:
                lower_factors.append(i)
        i += 1
    yield from reversed(lower_factors)


def _define_rle_allele(
    allele, length, repeat_subunit_length, rle_seq_limit, extended_alt_seq
):
    # Otherwise, create the Allele as an RLE
    allele.state = models.ReferenceLengthExpression(
        length=length, repeatSubunitLength=repeat_subunit_length
    )

    if (rle_seq_limit and length <= rle_seq_limit) or (rle_seq_limit is None):
        allele.state.sequence = models.sequenceString(extended_alt_seq)

    return allele


def _is_valid_cycle(template_start, template, target):
    cycle = itertools.cycle(template[template_start:])
    for char in target[len(template) :]:  # noqa: SIM110
        if char != next(cycle):
            return False
    return True


# TODO _normalize_genotype?


def _normalize_cis_phased_block(
    o,
    data_proxy: _DataProxy | None = None,  # noqa: ARG001
):
    o.members = sorted(o.members, key=ga4gh_digest)
    return o


handlers = {
    "Allele": _normalize_allele,
    "CisPhasedBlock": _normalize_cis_phased_block,
}


def normalize(vo, data_proxy: _DataProxy | None = None, **kwargs):
    """Normalize given vrs object, regardless of type

    :param vo:
    :param data_proxy: GA4GH sequence dataproxy instance, if needed
    :keyword rle_seq_limit: If RLE is set as the new state, set the limit for the length
        of the `sequence`. To exclude `state.sequence`, set to 0.
    :return: normalized object, or unmodified input object if the normalization algorithm
        does not provide normalization steps for the given type.
    :raise TypeError: if given object isn't a pydantic.BaseModel
    """
    if not is_pydantic_instance(vo):
        msg = f"Object is of class {vo.__class__} (with parents {vo.__class__.__mro__}), but normalize requires objects that inherit from pydantic.BaseModel."
        raise TypeError(msg)
    vo_type = vo.type

    if vo_type in handlers:
        handler = handlers[vo_type]
        return handler(vo, data_proxy, **kwargs)

    # No handler for vo_type; pass-through unchanged
    return vo
