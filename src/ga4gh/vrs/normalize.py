"""Normalize VRS objects.

Allele normalization validates bounds, trims shared sequence, expands repeat
ambiguity, and chooses a literal or reference-length expression. Substitutions
and identity alleles return before expansion.

For reference ``CAAAG``, inserting ``A`` inside the run expands as follows::

    input:    [2, 2), literal "A"
    expanded: [1, 4), reference "AAA", alternate "AAAA"
    output:   ReferenceLengthExpression(length=4, repeatSubunitLength=1)

A substitution needs only trimming: ``ACGT -> ATGT`` at [0, 4) becomes
``C -> T`` at [1, 2), with a LiteralSequenceExpression.

Given an inline literal Allele and a data proxy serving its reference::

    from ga4gh.vrs import normalize

    compact = normalize(allele, data_proxy, rle_seq_limit=0)
    with_sequence = normalize(allele, data_proxy, rle_seq_limit=None)

0 omits the optional RLE sequence; None retains it without a length limit.
A positive limit retains the complete sequence only if it fits. Use the original
literal input for each call: existing RLE states are returned unchanged.

See https://vrs.ga4gh.org/en/stable/impl-guide/normalization.html
"""

import itertools
import logging
from collections.abc import Iterator
from enum import IntEnum
from typing import NamedTuple

from bioutils.normalize import NormalizationMode
from bioutils.normalize import normalize as _normalize
from pydantic.main import BaseModel

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


def _normalize_allele(
    input_allele: models.Allele, data_proxy: _DataProxy, rle_seq_limit: int | None = 50
) -> models.Allele:
    """Normalize an Allele using fully justified normalization.

    Adapted from NCBI's VOCA, this expands ambiguous representations over the entire
    region of ambiguity so they can be compared with other alleles.

    This function assumes that IRIs are dereferenced, providing a `SequenceReference` as
    the `allele.location.sequenceReference`. If a `SequenceReference` is not provided,
    the allele will be returned as is with no normalization.

    Does not attempt to normalize Alleles with definite ranges or
    non-LiteralSequenceExpression states and will instead return the `input_allele`.

    If attempting to re-normalize a ReferenceLengthExpression allele, use
    denormalize_reference_length_expression to construct a LiteralSequenceExpression
    allele and then call normalize() on that LiteralSequenceExpression allele.

    See:
    https://vrs.ga4gh.org/en/2.0/conventions/normalization.html#literalsequenceexpression-alleles

    :param input_allele: Input VRS Allele object
    :param data_proxy: SeqRepo dataproxy
    :param rle_seq_limit: If RLE is set as the new state, set the limit for the length
        of the `sequence`.
        To exclude `sequence` from the response, set to 0.
        For no limit, set to `None`.
    :return: Normalized allele, or the input if normalization is unsupported.
    :raises DataProxyValidationError: If the allele location is out of bounds on its
        sequence
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

    # Reject locations that do not exist on the sequence before anything is fetched,
    # since out-of-range fetches may be silently truncated by the sequence backend.
    # Done before the early returns below, which skip definite ranges.
    data_proxy.validate_location_bounds(
        alias, input_allele.location.start, input_allele.location.end
    )

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


def _trim_for_normalization(
    ref_seq: SequenceProxy,
    ival: tuple[int, int],
    alleles: tuple[None, str],
    start: LocationPos,
    end: LocationPos,
) -> tuple[tuple[int, int], tuple[str | None, str]]:
    """Trim common prefix and suffix from the intervals.

    Return the trimmed interval and trimmed alleles:
        ((trim_start, trim_end), (trim_ref, trim_alt))

    If the alleles are the same, return the original interval and alleles. The first
    allele (ref) will be populated by bioutils.

    :param ref_seq: Reference sequence proxy.
    :param ival: Input start and end coordinates.
    :param alleles: Reference placeholder and alternate sequence.
    :param start: Representative start position.
    :param end: Representative end position.
    :return: Trimmed interval and allele sequences.
    :raises ValueError: If trimming fails for a non-reference allele.
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


def _set_location_from_interval(
    allele: models.Allele,
    ival: tuple[int, int],
    start_pos_type: PosType,
    end_pos_type: PosType,
) -> None:
    """Update ``allele`` start and end location

    :param allele: Allele to modify in place.
    :param ival: Start and end coordinates.
    :param start_pos_type: Representation of the original start position.
    :param end_pos_type: Representation of the original end position.
    """
    allele.location.start = _get_new_allele_location_pos(ival[0], start_pos_type)
    allele.location.end = _get_new_allele_location_pos(ival[1], end_pos_type)


def denormalize_reference_length_expression(
    ref_seq: str,
    repeat_subunit_length: int,
    alt_length: int,
) -> str:
    """Expand a ReferenceLengthExpression into its alternate sequence.

    Repeat the reference subunit to alt_length characters, including a trailing partial
    copy when needed.

    e.g. "ACGT" (repeat_subunit_length=4, length=10) -> "ACGTACGTAC" (ACGT ACGT AC)

    :param ref_seq: The reference sequence to be denormalized.
    :param repeat_subunit_length: The length of the repeat subunit.
    :param alt_length: The length of the alternate sequence that was compacted during normalization.

    :return: Expanded alternate sequence.
    :raises ValueError: If the subunit length is negative or exceeds the reference length.
    :raises ZeroDivisionError: If the repeat subunit length is zero.
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


def _factor_gen(n: int) -> Iterator[int]:
    """Yield all factors of an integer in ascending order.

    The first factor that passes `_is_valid_cycle` becomes the ``repeatSubunitLength``,
    so ascending order selects the smallest valid factor (VRS 2.1 normalization). VRS
    2.0 selected the largest; see ``PrevVrsVersion.V2_0`` for computing 2.0 digests.

    :param n: Positive integer to factor. Nonpositive inputs yield no values.
    :return: Iterator yielding factors from smallest to largest.
    """
    upper_factors = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            yield i
            if n // i != i:
                upper_factors.append(n // i)
        i += 1
    yield from reversed(upper_factors)


def _define_rle_allele(
    allele: models.Allele,
    length: int,
    repeat_subunit_length: int,
    rle_seq_limit: int | None,
    extended_alt_seq: str,
) -> models.Allele:
    """Replace an allele state with a reference length expression.

    :param allele: Allele to modify in place.
    :param length: Alternate sequence length.
    :param repeat_subunit_length: Length of the repeating unit.
    :param rle_seq_limit: Maximum sequence length to retain, or None for no limit.
    :param extended_alt_seq: Expanded alternate sequence.
    :return: Updated allele.
    """
    allele.state = models.ReferenceLengthExpression(
        length=length, repeatSubunitLength=repeat_subunit_length
    )

    if (rle_seq_limit and length <= rle_seq_limit) or (rle_seq_limit is None):
        allele.state.sequence = models.sequenceString(extended_alt_seq)

    return allele


def _is_valid_cycle(template_start: int, template: str, target: str) -> bool:
    """Check whether the target suffix repeats the template cycle.

    :param template_start: Start of a nonempty repeating unit within the template.
    :param template: Reference sequence containing the repeating unit.
    :param target: Alternate sequence whose suffix is checked.
    :return: Whether all suffix characters match the cycle. The template-length
        prefix is not compared; an empty suffix matches.
    :raises StopIteration: If the cycle is empty and the target suffix is not.
    """
    cycle = itertools.cycle(template[template_start:])
    for char in target[len(template) :]:  # noqa: SIM110
        if char != next(cycle):
            return False
    return True


# TODO _normalize_genotype?


def _normalize_cis_phased_block(
    o: models.CisPhasedBlock,
    data_proxy: _DataProxy | None = None,  # noqa: ARG001
) -> models.CisPhasedBlock:
    """Sort block members by their GA4GH digests in place.

    :param o: Block to normalize.
    :param data_proxy: Unused sequence data proxy.
    :return: Updated block.
    """
    o.members = sorted(o.members, key=ga4gh_digest)
    return o


handlers = {
    "Allele": _normalize_allele,
    "CisPhasedBlock": _normalize_cis_phased_block,
}


def normalize(
    vo: BaseModel, data_proxy: _DataProxy | None = None, **kwargs: int | None
) -> BaseModel:
    """Normalize given vrs object, regardless of type

    :param vo: VRS object to normalize.
    :param kwargs: Allele handler options. rle_seq_limit limits retained RLE sequence
        length; 0 omits it and None removes the limit. CisPhasedBlock accepts
        no extra options.
    :param data_proxy: Required for supported Allele normalization; unused by
        CisPhasedBlock normalization.
    :return: Normalized object of the input type. Alleles may be copied;
        CisPhasedBlock members are sorted in place. Unsupported types are
        returned unchanged. The input must expose a VRS type attribute.
    :raises TypeError: if given object isn't a pydantic.BaseModel
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
