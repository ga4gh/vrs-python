"""Convert HGVS expressions and coordinates.

HGVS uses 1-based positions; VRS uses interbase bounds. Omitting the accession::

    HGVS                            VRS start        VRS end
    g.100A>T                        99               100
    g.100_101insT                    100              100
    g.(100_110)_(200_210)del          Range([99,109])  Range([200,210])

Insertions use the boundary after the start residue for both bounds. For ranges,
starts shift by one and ends stay unchanged. AlleleTranslator rejects uncertain
ranges; use CnvTranslator for uncertain del/dup copy-number expressions.

With a configured data_proxy and accessible UTA database::

    from ga4gh.vrs.extras.translator import AlleleTranslator

    translator = AlleleTranslator(data_proxy)
    try:
        allele = translator.from_hgvs("NC_000019.10:g.44908822C>T")
        expressions = translator.to_hgvs(allele, namespace="refseq")
    finally:
        translator.hgvs_tools.close()  # Closes UTA, not the supplied data proxy.

Coding coordinates pass through c_to_n on import and n_to_c on export.
"""

import logging
import re
from types import TracebackType
from typing import Literal, Self, TypeGuard

import hgvs
import hgvs.dataproviders.uta
import hgvs.edit
import hgvs.location
import hgvs.normalizer
import hgvs.parser
import hgvs.posedit
import hgvs.variantmapper
from hgvs.sequencevariant import SequenceVariant as HgvsSequenceVariant

from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import _DataProxy

_logger = logging.getLogger(__name__)

# HGVS uses 1-based inclusive coordinates; VRS uses 0-based interbase. For the
# outer-left side of a variant interval, HGVS position N corresponds to VRS
# position N-1; for the outer-right side, HGVS N corresponds to VRS N.
_Side = Literal["start", "end"]

_HgvsPos = (
    hgvs.location.SimplePosition
    | hgvs.location.BaseOffsetPosition
    | hgvs.location.Interval
)
_VrsPos = int | models.Range | None


def _is_uncertain_range(pos: _HgvsPos) -> TypeGuard[hgvs.location.Interval]:
    """Return True if ``pos`` is a nested :class:`hgvs.location.Interval` representing an uncertain range (e.g. ``(A_B)`` from ``(A_B)_(C_D)del``).

    Note that hgvs 2.0.0a0 also emits a non-uncertain ``Interval`` (with
    ``start == end``) for the exact side of a mixed expression like
    ``g.100_(200_?)del``, so the ``.uncertain`` check is required to distinguish
    real uncertain ranges from that wrapper shape.

    :param pos: HGVS position to inspect or convert.
    :return: Whether the position is an uncertain interval.
    """
    return isinstance(pos, hgvs.location.Interval) and pos.uncertain


def _has_intron_offset(pos: _HgvsPos) -> bool:
    """Return True if ``pos`` is a :class:`hgvs.location.BaseOffsetPosition` with a non-zero intron offset (e.g. ``100+5`` or ``200-3``).

    VRS :class:`SequenceLocation` coordinates are transcript-relative; a
    non-zero offset means the position lies inside an intron, which has no
    faithful single-integer representation in VRS. This helper lets the
    low-level converters refuse such inputs defensively; upstream callers
    typically reject at the variant level via :meth:`HgvsTools.is_intronic`.

    :param pos: HGVS position to inspect or convert.
    :return: Whether the position has a nonzero intron offset.
    """
    return isinstance(pos, hgvs.location.BaseOffsetPosition) and bool(pos.offset)


def _shift_hgvs_to_vrs(base: int | None, side: _Side) -> int | None:
    """Convert an HGVS 1-based inclusive coordinate to a VRS 0-based interbase coordinate, accounting for which side of the outer interval it is on.

    :param base: HGVS coordinate, or None if unknown.
    :param side: Start or end side of the variant interval.
    :return: Interbase coordinate, or None if unknown.
    """
    if base is None:
        return None
    return base - 1 if side == "start" else base


def _shift_vrs_to_hgvs(value: int | None, side: _Side) -> int | None:
    """Inverse of :func:`_shift_hgvs_to_vrs`.

    :param value: Value to convert.
    :param side: Start or end side of the variant interval.
    :return: HGVS coordinate, or None if unknown.
    """
    if value is None:
        return None
    return value + 1 if side == "start" else value


def _hgvs_pos_to_vrs(pos: _HgvsPos, side: _Side) -> _VrsPos:
    """Convert the ``pos.start`` or ``pos.end`` of a parsed hgvs variant to the corresponding VRS value for :attr:`models.SequenceLocation.start` / :attr:`models.SequenceLocation.end`.

    Handles three shapes of ``pos``:

    * :class:`hgvs.location.SimplePosition` / :class:`hgvs.location.BaseOffsetPosition`
      represents a certain position; returns an ``int``.
    * :class:`hgvs.location.Interval` with ``uncertain=True``, a nested
      uncertain range (e.g. ``(A_B)`` from ``(A_B)_(C_D)dup``); returns a
      :class:`models.Range`.
    * :class:`hgvs.location.Interval` with ``uncertain=False``, the
      ``start==end`` wrapper hgvs emits for the exact side of a mixed
      expression like ``g.100_(200_?)del``; unwrapped to an ``int``.

    Does **not** accept :class:`hgvs.location.BaseOffsetInterval` or any
    intronic :class:`~hgvs.location.BaseOffsetPosition` (non-zero ``offset``):
    those have no faithful single-integer VRS representation, and this helper
    raises rather than silently dropping offset information. Callers that want
    a variant-level error (rather than a position-level one) may reject with
    :meth:`HgvsTools.is_intronic` first.

    :param side: ``"start"`` if ``pos`` is the outer-left side of the variant
        (HGVS→VRS subtracts 1), ``"end"`` for the outer-right side (no shift).
    :return: An ``int`` for a certain position, a :class:`models.Range` for an
        uncertain range, or ``None`` if the position is fully unknown.
    :raises TypeError: if ``pos`` is a :class:`~hgvs.location.BaseOffsetInterval`
        (c./n. transcript range), which has no VRS equivalent.
    :raises ValueError: if ``pos`` is an intronic
        :class:`~hgvs.location.BaseOffsetPosition` (non-zero ``offset``), or a
        non-uncertain interval has different endpoints.

    :param pos: HGVS position to inspect or convert.
    """
    if isinstance(pos, hgvs.location.Interval):
        if isinstance(pos, hgvs.location.BaseOffsetInterval):
            msg = (
                "BaseOffsetInterval (c./n. transcript range) is not "
                "representable as a VRS SequenceLocation coordinate"
            )
            raise TypeError(msg)
        if pos.uncertain:
            lo = _shift_hgvs_to_vrs(pos.start.base, side)
            hi = _shift_hgvs_to_vrs(pos.end.base, side)
            return models.Range(root=[lo, hi])
        # Exact side of a mixed expression like g.100_(200_?)del: hgvs wraps
        # the certain endpoint in a non-uncertain Interval with start==end so
        # both outer sides share a type. Unwrap to the underlying base. The
        # inner positions are plain SimplePositions here (this wrapper shape
        # only appears for g. expressions), so no intron-offset check needed.
        if pos.start.base != pos.end.base:
            msg = (
                f"Unexpected non-uncertain Interval with start.base "
                f"({pos.start.base}) != end.base ({pos.end.base}); only the "
                f"hgvs mixed-endpoint wrapper shape is supported here"
            )
            raise ValueError(msg)
        return _shift_hgvs_to_vrs(pos.start.base, side)
    if _has_intron_offset(pos):
        msg = (
            "Intronic position (non-zero BaseOffsetPosition.offset) is not "
            "representable as a VRS SequenceLocation coordinate"
        )
        raise ValueError(msg)
    return _shift_hgvs_to_vrs(pos.base, side)


def _vrs_pos_to_hgvs(
    value: int | models.Range | list[int | None] | None,
    side: _Side,
    *,
    as_interval: bool = False,
) -> hgvs.location.SimplePosition | hgvs.location.Interval:
    """Convert a VRS :attr:`models.SequenceLocation.start` / ``end`` value to the corresponding hgvs position object.

    For an ``int``, returns a :class:`hgvs.location.SimplePosition`. For a
    :class:`models.Range` (or equivalent 2-element list), returns a nested
    :class:`hgvs.location.Interval` with ``uncertain=True``. ``None`` bounds
    within a Range become ``SimplePosition(base=None)`` (rendered as ``?``).

    When ``as_interval`` is True, an ``int`` value is wrapped in a non-uncertain
    ``Interval`` with ``start == end``. This matches the hgvs parser's
    representation of mixed certain/uncertain outer intervals (e.g.
    ``g.N_(M_?)del``), where both sides of the outer interval must have the
    same kind of inner position for the hgvs formatter to work.

    :param value: Value to convert.
    :param side: Start or end side of the variant interval.
    :param as_interval: Wrap an exact position in an interval when true.
    :return: HGVS position or interval.
    """
    if isinstance(value, models.Range):
        value = value.root
    if isinstance(value, list):
        lo = _shift_vrs_to_hgvs(value[0], side)
        hi = _shift_vrs_to_hgvs(value[1], side)
        return hgvs.location.Interval(
            start=hgvs.location.SimplePosition(base=lo),
            end=hgvs.location.SimplePosition(base=hi),
            uncertain=True,
        )
    base = _shift_vrs_to_hgvs(value, side)
    if as_interval:
        return hgvs.location.Interval(
            start=hgvs.location.SimplePosition(base=base),
            end=hgvs.location.SimplePosition(base=base),
        )
    return hgvs.location.SimplePosition(base=base)


class HgvsTools:
    """Provide tools for working with HGVS notation.

    Use a context manager or call ``close()`` to close this instance's UTA
    connection. The caller owns the supplied data proxy, which stays open.

    Attributes:
        parser (hgvs.parser.Parser): The HGVS parser object.
        uta_conn: The UTA (Universal Transcript Archive) connection object.
        normalizer: The HGVS normalizer object.
        variant_mapper: The HGVS variant mapper object.
        data_proxy (_DataProxy): The data proxy object.

    """

    hgvs_re = re.compile(r"[^:]+:[cgmnpr]\.")

    def __init__(self, data_proxy: _DataProxy | None = None) -> None:
        """Create the HGVS parser, normalizer, variant mapper, and UTA connection.

        Close the UTA connection if normalizer or variant mapper setup fails.
        If cleanup also fails, raise its error with the setup error chained.

        :param data_proxy: Optional caller-owned GA4GH data proxy.
        """
        self.parser = hgvs.parser.Parser()
        self.uta_conn = hgvs.dataproviders.uta.connect()
        self._closed = False
        try:
            self.normalizer = hgvs.normalizer.Normalizer(self.uta_conn, validate=True)
            self.variant_mapper = hgvs.variantmapper.VariantMapper(self.uta_conn)
        except BaseException:
            self.close()
            raise
        self.data_proxy = data_proxy

    def close(self) -> None:
        """Close this instance's UTA connection.

        Repeated calls after success do nothing. Cleanup errors propagate;
        callers can retry. Keep all resource references and leave the supplied
        data proxy open. Operations that need UTA are unsupported after closure.
        """
        if not self._closed:
            self.uta_conn.close()
            self._closed = True

    def __enter__(self) -> Self:
        """Return this instance for use in a context manager.

        :return: This instance.
        :raises RuntimeError: If already closed.
        """
        if self._closed:
            msg = "HgvsTools is closed"
            raise RuntimeError(msg)
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Close the owned UTA connection without suppressing exceptions.

        If cleanup fails, raise its error with any context body error chained.

        :param exc_type: Exception type from the context body, if any.
        :param exc_value: Exception from the context body, if any.
        :param traceback: Traceback from the context body, if any.
        """
        self.close()

    # convenience methods for hgvs parsing, normalization, and some mappings
    def parse(self, hgvs_str: str) -> HgvsSequenceVariant | None:
        """Parse the given HGVS string and returns the corresponding variant.

        :param hgvs_str: HGVS expression to parse.
        :return: Parsed variant, or None if the expression prefix is unsupported.
        :raises hgvs.exceptions.HGVSParseError: If a supported expression has invalid syntax.
        """
        if not self.hgvs_re.match(hgvs_str):
            return None
        return self.parser.parse_hgvs_variant(hgvs_str)

    def has_uncertain_range(self, sv: HgvsSequenceVariant) -> bool:
        """Return True if either endpoint of ``sv.posedit.pos`` is a nested uncertain :class:`hgvs.location.Interval` (e.g. the ``(A_B)`` sides of ``(A_B)_(C_D)del``). Surfaces the module-private :func:`_is_uncertain_range` check to external callers without requiring them to import an underscore helper across module boundaries.

        :param sv: Parsed HGVS sequence variant.
        :return: Whether either endpoint is an uncertain range.
        """
        return _is_uncertain_range(sv.posedit.pos.start) or _is_uncertain_range(
            sv.posedit.pos.end
        )

    def is_intronic(self, sv: HgvsSequenceVariant) -> bool:
        """Check if the given SequenceVariant is intronic.

        :param sv: Parsed HGVS sequence variant.
        :return: Whether either endpoint is intronic.
        """
        # Coding (c.) intervals are BaseOffsetInterval, but RNA (r.) intervals
        # are plain Interval even when their endpoints carry a base-offset
        # position, so inspect the endpoints directly rather than the interval
        # type. Non-offset coordinates (g./p.) use position types without
        # is_intronic and are correctly reported as not intronic.
        pos = sv.posedit.pos
        start = getattr(pos, "start", None)
        end = getattr(pos, "end", None)
        return (
            isinstance(start, hgvs.location.BaseOffsetPosition) and start.is_intronic
        ) or (isinstance(end, hgvs.location.BaseOffsetPosition) and end.is_intronic)

    def get_edit_type(self, sv: HgvsSequenceVariant | None) -> str | None:
        """Safely extract ``type`` property from SequenceVariant

        :param sv: Parsed HGVS sequence variant.
        :return: Edit type, or None if no edit is present.
        """
        if sv is None or sv.posedit is None or sv.posedit.edit is None:
            return None
        return sv.posedit.edit.type

    def get_position_and_state(
        self, sv: HgvsSequenceVariant
    ) -> tuple[_VrsPos, _VrsPos, str]:
        """Get the details of a sequence variant.

        :param sv: Parsed HGVS sequence variant.
        :return: VRS start, end, and literal alternate sequence.
        :raises ValueError: If the edit type or an intronic position is unsupported.
        :raises TypeError: If a transcript range cannot be represented in VRS.
        :raises AttributeError: If an operation needing sequence data is called without a data proxy.
        """
        pos = sv.posedit.pos
        edit_type = sv.posedit.edit.type

        # hgvs insertion semantics: the insertion sits between pos.start and
        # pos.start+1, so both VRS bounds take pos.start.base with no -1
        # shift. Doesn't fit the side="start"/"end" convention, so bypass
        # _hgvs_pos_to_vrs here and read the attribute directly.
        if edit_type == "ins":
            start = end = pos.start.base
            state = sv.posedit.edit.alt
            return start, end, state

        # For every other edit type, pos.start and pos.end convert to VRS
        # coordinates the same way. _hgvs_pos_to_vrs raises typed errors for
        # intronic offsets and BaseOffsetInterval inputs, so the arithmetic
        # below is safe against those shapes without explicit guards here.
        start = _hgvs_pos_to_vrs(pos.start, side="start")
        end = _hgvs_pos_to_vrs(pos.end, side="end")

        if edit_type in ("sub", "del", "delins", "identity"):
            if edit_type == "identity":
                # Checked before fetching, which a backend may silently truncate
                self.data_proxy.validate_location_bounds(sv.ac, start, end)
                state = self.data_proxy.get_sequence(sv.ac, start=start, end=end)
            else:
                state = sv.posedit.edit.alt or ""

        elif edit_type == "dup":
            self.data_proxy.validate_location_bounds(sv.ac, start, end)
            ref = self.data_proxy.get_sequence(sv.ac, start=start, end=end)
            state = ref + ref

        else:
            msg = f"HGVS variant type {edit_type} is unsupported"
            raise ValueError(msg)

        return start, end, state

    def extract_allele_values(self, hgvs_expr: str) -> dict | None:
        """Extract allele construction fields from an HGVS expression.

        :param hgvs_expr: HGVS expression to translate.
        :return: Allele fields, or None for an unsupported expression prefix or
            unresolved accession. Invalid HGVS syntax raises HGVSParseError.
        :raises ValueError: If the variant is intronic or its edit type is unsupported.
        :raises AttributeError: If an operation needing sequence data is called without a data proxy.
        :raises TypeError: If a transcript range cannot be represented in VRS.
        :raises hgvs.exceptions.HGVSParseError: If supported HGVS syntax is malformed.
        """
        sv = self.parse(hgvs_expr)
        if not sv:
            return None

        if self.is_intronic(sv):
            msg = "Intronic HGVS variants are not supported"
            raise ValueError(msg)

        refget_accession = self.data_proxy.derive_refget_accession(sv.ac)
        if not refget_accession:
            return None

        # translate coding coordinates to positional coordinates, if necessary
        if sv.type == "c":
            sv = self.c_to_n(sv)

        (start, end, state) = self.get_position_and_state(sv)

        return {
            "sequence_id": sv.ac,
            "refget_accession": refget_accession,
            "start": start,
            "end": end,
            "literal_sequence": state,
        }

    def build_cnv_location(
        self, sv: HgvsSequenceVariant, refget_accession: str
    ) -> models.SequenceLocation:
        """Build a VRS :class:`SequenceLocation` from a parsed hgvs CNV variant.

        CDS-relative (``c.``) inputs are translated to transcript-relative
        (``n.``) coordinates first; without that step the VRS location would
        point into the 5' UTR. Handles both certain positions and uncertain-
        range endpoints; the resulting ``start`` / ``end`` are ints for exact
        positions and :class:`models.Range` values for uncertain ranges.

        :param sv: Parsed HGVS sequence variant.
        :param refget_accession: Refget accession for the location sequence.
        :return: Sequence location with translated coordinates.
        :raises TypeError: If a transcript range cannot be represented in VRS.
        :raises ValueError: If a position is intronic or an exact interval has unequal endpoints.
        """
        if sv.type == "c":
            sv = self.c_to_n(sv)
        return models.SequenceLocation(
            sequenceReference=models.SequenceReference(
                refgetAccession=refget_accession
            ),
            start=_hgvs_pos_to_vrs(sv.posedit.pos.start, side="start"),
            end=_hgvs_pos_to_vrs(sv.posedit.pos.end, side="end"),
        )

    def from_allele(
        self, vo: models.Allele | None, namespace: str | None = None
    ) -> list[str]:
        """Generate a *list* of HGVS expressions for VRS Allele.

        If `namespace` is not None, returns HGVS strings for the
        specified namespace.

        If `namespace` is None, returns HGVS strings for all alias
        translations.

        If no alias translations are available, an empty list is
        returned.

        If the VRS object cannot be expressed as HGVS, raises ValueError.

        :param vo: VRS object to translate.
        :param namespace: Optional namespace used to filter sequence aliases.
        :return: Distinct HGVS expressions for supported matching aliases, or an empty
            list for None. Protein accessions raise ValueError.
        :raises ValueError: If the allele or sequence reference cannot be represented.
        :raises AttributeError: If an operation needing sequence data is called without a data proxy.
        """
        if vo is None:
            return []
        if not isinstance(vo, models.Allele):
            msg = "VRS object must be an Allele"
            raise ValueError(msg)  # noqa: TRY004
        if vo.location is None:
            msg = "VRS allele must have a location"
            raise ValueError(msg)

        refget_accession = vo.location.get_refget_accession()
        if refget_accession is None:
            msg = "VRS allele location must have a sequence reference"
            raise ValueError(msg)

        sequence = f"ga4gh:{refget_accession}"
        aliases = self.data_proxy.translate_sequence_identifier(sequence, namespace)

        hgvs_exprs = []
        for alias in aliases:
            ns, accession = alias.split(":")
            # skip GRCh accessions unless specifically requested
            # because they are ambiguous without their namespace,
            # which can't be included in HGVS expressions
            # TODO: use default_assembly_name here
            if ns.startswith("GRC") and namespace is None:
                continue

            if not (
                any(
                    accession.startswith(pfx)
                    for pfx in (
                        "NM",
                        "NP",
                        "NC",
                        "NG",
                        "NR",
                        "NW",
                        "NT",
                        "XM",
                        "XR",
                        "XP",
                    )
                )
            ):
                continue

            sequence_type = self.data_proxy.extract_sequence_type(alias)

            # create the hgvs expression object
            var = self._to_sequence_variant(vo, sequence_type, sequence, accession)
            hgvs_exprs += [str(var)]

        return list(set(hgvs_exprs))

    def _to_sequence_variant(
        self,
        vo: models.Allele,
        sequence_type: str | None,
        sequence: str,
        accession: str,
    ) -> HgvsSequenceVariant:
        """Create a SequenceVariant object from an Allele object.

        :param vo: VRS object to translate.
        :param sequence_type: HGVS sequence type, such as g, n, or c; None when the
            alias prefix is unknown.
        :param sequence: Sequence identifier used to retrieve reference bases.
        :param accession: Accession for the HGVS expression.
        :return: HGVS variant representing the allele. Missing UTA data is logged
            and leaves the constructed variant unnormalized.
        :raises ValueError: If the sequence type is protein.
        """
        # build interval and edit depending on sequence type
        if sequence_type == "p":
            msg = "Only nucleic acid variation is currently supported"
            raise ValueError(msg)
            # ival = hgvs.location.Interval(start=start, end=end)
            # edit = hgvs.edit.AARefAlt(ref=None, alt=vo.state.sequence)
        start, end = vo.location.start, vo.location.end
        # ib: 0 1 2 3 4 5
        #  h:  1 2 3 4 5
        if isinstance(start, models.Range) or isinstance(end, models.Range):
            # Uncertain-range bounds: emit nested uncertain Intervals. The
            # reference sequence is unknown; pass "" so hgvs formats a bare
            # "del"/"delins"-style edit rather than erroring on ref=None. If
            # only one side is uncertain, wrap the certain side in a non-
            # uncertain Interval(start==end) so hgvs's outer Interval format
            # can compare start and end (it asserts matching types).
            ref = ""
            ival = hgvs.location.Interval(
                start=_vrs_pos_to_hgvs(start, side="start", as_interval=True),
                end=_vrs_pos_to_hgvs(end, side="end", as_interval=True),
            )
        else:
            if start == end:  # insert: hgvs uses *exclusive coords*
                ref = None
                end += 1
            else:  # else: hgvs uses *inclusive coords*
                ref = self.data_proxy.get_sequence(sequence, start, end)
                start += 1

            ival = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(start),
                end=hgvs.location.SimplePosition(end),
            )
        alt = str(vo.state.sequence.root) or None  # "" => None
        edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)

        posedit = hgvs.posedit.PosEdit(pos=ival, edit=edit)
        var = HgvsSequenceVariant(
            ac=accession,
            # at this point, use `n.` because the positions are absolute (not CDS),
            # this will subsequently be converted back to `c.` after hgvs normalization
            type="n" if sequence_type == "c" else sequence_type,
            posedit=posedit,
        )

        try:
            # if the namespace is GRC, can't normalize, since hgvs can't deal with it
            parsed = self.parse(str(var))
            var = self.normalize(parsed)

            # if sequence_type is coding, convert from "n." to "c." before continuing
            if sequence_type == "c":
                var = self.n_to_c(var)

        except hgvs.exceptions.HGVSDataNotAvailableError:
            _logger.warning("No data found for accession %s", accession)

        return var

    def normalize(self, hgvs: HgvsSequenceVariant) -> HgvsSequenceVariant:
        """Perform hgvs library `normalize()` method

        :param hgvs: HGVS variant to normalize or map.
        :return: Normalized HGVS variant.
        :raises hgvs.exceptions.HGVSError: If the HGVS operation fails validation or required data is unavailable.
        """
        return self.normalizer.normalize(hgvs)

    def n_to_c(self, hgvs: HgvsSequenceVariant) -> HgvsSequenceVariant:
        """Perform hgvs library `n_to_c` method

        :param hgvs: HGVS variant to normalize or map.
        :return: Variant in coding coordinates.
        :raises hgvs.exceptions.HGVSError: If the HGVS operation fails validation or required data is unavailable.
        """
        return self.variant_mapper.n_to_c(hgvs)

    def c_to_n(self, hgvs: HgvsSequenceVariant) -> HgvsSequenceVariant:
        """Perform hgvs library `c_to_n` method

        :param hgvs: HGVS variant to normalize or map.
        :return: Variant in transcript coordinates.
        :raises hgvs.exceptions.HGVSError: If the HGVS operation fails validation or required data is unavailable.
        """
        return self.variant_mapper.c_to_n(hgvs)


if __name__ == "__main__":
    import os

    from ga4gh.vrs.dataproxy import create_dataproxy

    seqrepo_uri = os.environ.get(
        "SEQREPO_URI", "seqrepo+file:///usr/local/share/seqrepo/latest"
    )
    if seqrepo_uri is None:
        msg = "SEQREPO_URI environment variable must be set to a valid seqrepo URI"
        raise ValueError(msg)

    allele_dict = {
        "id": "ga4gh:VA.SmhjExyBS8GxicngJxQLJ8Ww5GrBQk40",
        "type": "Allele",
        "digest": "SmhjExyBS8GxicngJxQLJ8Ww5GrBQk40",
        "location": {
            "id": "ga4gh:SL.HyEEyTvdJrfsElsTXRG-43Y_tID4QYFt",
            "type": "SequenceLocation",
            "digest": "HyEEyTvdJrfsElsTXRG-43Y_tID4QYFt",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.gg9z_DUcgwhfRz29BS3zq2jim59tOrA9",
            },
            "start": 874,
            "end": 875,
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "G"},
    }

    vrs_allele = models.Allele.model_validate(allele_dict)
    dp = create_dataproxy(seqrepo_uri)
    with HgvsTools(dp) as hgvs_tools:
        hgvs_expr = hgvs_tools.from_allele(vrs_allele, namespace="refseq")

    print(hgvs_expr)  # noqa: T201
