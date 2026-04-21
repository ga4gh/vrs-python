"""Provide extra tools for working with hgvs expressions."""

import logging
import re
from typing import Literal, TypeGuard

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


def _is_uncertain_range(pos: _HgvsPos) -> TypeGuard[hgvs.location.Interval]:
    """Return True if ``pos`` is a nested :class:`hgvs.location.Interval`
    representing an uncertain range (e.g. ``(A_B)`` from ``(A_B)_(C_D)del``).

    Note that hgvs 2.0.0a0 also emits a non-uncertain ``Interval`` (with
    ``start == end``) for the exact side of a mixed expression like
    ``g.100_(200_?)del``, so the ``.uncertain`` check is required to distinguish
    real uncertain ranges from that wrapper shape.
    """
    return isinstance(pos, hgvs.location.Interval) and pos.uncertain


def _has_intron_offset(pos: _HgvsPos) -> bool:
    """Return True if ``pos`` is a :class:`hgvs.location.BaseOffsetPosition`
    with a non-zero intron offset (e.g. ``100+5`` or ``200-3``).

    VRS :class:`SequenceLocation` coordinates are transcript-relative; a
    non-zero offset means the position lies inside an intron, which has no
    faithful single-integer representation in VRS. This helper lets the
    low-level converters refuse such inputs defensively; upstream callers
    typically reject at the variant level via :meth:`HgvsTools.is_intronic`.
    """
    return isinstance(pos, hgvs.location.BaseOffsetPosition) and bool(pos.offset)


def _shift_hgvs_to_vrs(base: int | None, side: _Side) -> int | None:
    """Convert an HGVS 1-based inclusive coordinate to a VRS 0-based interbase
    coordinate, accounting for which side of the outer interval it is on.
    """
    if base is None:
        return None
    return base - 1 if side == "start" else base


def _shift_vrs_to_hgvs(value: int | None, side: _Side) -> int | None:
    """Inverse of :func:`_shift_hgvs_to_vrs`."""
    if value is None:
        return None
    return value + 1 if side == "start" else value


def _hgvs_pos_to_vrs(pos: _HgvsPos, side: _Side) -> int | models.Range | None:
    """Convert the ``pos.start`` or ``pos.end`` of a parsed hgvs variant to the
    corresponding VRS value for :attr:`models.SequenceLocation.start` /
    :attr:`models.SequenceLocation.end`.

    Handles three shapes of ``pos``:

    * :class:`hgvs.location.SimplePosition` / :class:`hgvs.location.BaseOffsetPosition`
      — a certain position; returns an ``int``.
    * :class:`hgvs.location.Interval` with ``uncertain=True`` — a nested
      uncertain range (e.g. ``(A_B)`` from ``(A_B)_(C_D)dup``); returns a
      :class:`models.Range`.
    * :class:`hgvs.location.Interval` with ``uncertain=False`` — the
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
    :returns: An ``int`` for a certain position, a :class:`models.Range` for an
        uncertain range, or ``None`` if the position is fully unknown.
    :raises TypeError: if ``pos`` is a :class:`~hgvs.location.BaseOffsetInterval`
        (c./n. transcript range), which has no VRS equivalent.
    :raises ValueError: if ``pos`` is an intronic
        :class:`~hgvs.location.BaseOffsetPosition` (non-zero ``offset``).
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
    """Convert a VRS :attr:`models.SequenceLocation.start` / ``end`` value to
    the corresponding hgvs position object.

    For an ``int``, returns a :class:`hgvs.location.SimplePosition`. For a
    :class:`models.Range` (or equivalent 2-element list), returns a nested
    :class:`hgvs.location.Interval` with ``uncertain=True``. ``None`` bounds
    within a Range become ``SimplePosition(base=None)`` (rendered as ``?``).

    When ``as_interval`` is True, an ``int`` value is wrapped in a non-uncertain
    ``Interval`` with ``start == end``. This matches the hgvs parser's
    representation of mixed certain/uncertain outer intervals (e.g.
    ``g.N_(M_?)del``), where both sides of the outer interval must have the
    same kind of inner position for the hgvs formatter to work.
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
    """A convenience class that exposes only the tools needed by vrs-python for working with HGVS (Human Genome Variation Society) notation.

    Attributes:
        parser (hgvs.parser.Parser): The HGVS parser object.
        uta_conn: The UTA (Universal Transcript Archive) connection object.
        normalizer: The HGVS normalizer object.
        variant_mapper: The HGVS variant mapper object.
        data_proxy (_DataProxy): The data proxy object.

    """

    hgvs_re = re.compile(r"[^:]+:[cgmnpr]\.")

    def __init__(self, data_proxy: _DataProxy | None = None) -> None:
        """Initialize object.

        :param data_proxy: GA4GH data proxy instance
        """
        self.parser = hgvs.parser.Parser()
        self.uta_conn = hgvs.dataproviders.uta.connect()
        self.normalizer = hgvs.normalizer.Normalizer(self.uta_conn, validate=True)
        self.variant_mapper = hgvs.variantmapper.VariantMapper(self.uta_conn)
        self.data_proxy = data_proxy

    def close(self) -> None:
        """Close resources.

        # TODO These should only be closed if they are owned by this instance
        """
        self.normalizer = None
        self.variant_mapper = None
        self.data_proxy = None
        if self.uta_conn is not None:
            self.uta_conn.close()

    # convenience methods for hgvs parsing, normalization, and some mappings
    def parse(self, hgvs_str: str) -> HgvsSequenceVariant | None:
        """Parse the given HGVS string and returns the corresponding variant.

        Args:
            hgvs_str (str): The HGVS string to parse.

        Returns:
            Variant: The parsed variant object, or None if the HGVS string is invalid.

        """
        if not self.hgvs_re.match(hgvs_str):
            return None
        return self.parser.parse_hgvs_variant(hgvs_str)

    def has_uncertain_range(self, sv: HgvsSequenceVariant) -> bool:
        """Return True if either endpoint of ``sv.posedit.pos`` is a nested
        uncertain :class:`hgvs.location.Interval` (e.g. the ``(A_B)`` sides of
        ``(A_B)_(C_D)del``). Surfaces the module-private
        :func:`_is_uncertain_range` check to external callers without
        requiring them to import an underscore helper across module boundaries.
        """
        return _is_uncertain_range(sv.posedit.pos.start) or _is_uncertain_range(
            sv.posedit.pos.end
        )

    def is_intronic(self, sv: HgvsSequenceVariant) -> bool:
        """Check if the given SequenceVariant is intronic.

        Args:
            sv (hgvs.sequencevariant.SequenceVariant): The SequenceVariant to check.

        Returns:
            bool: True if the SequenceVariant is intronic, False otherwise.

        """
        if isinstance(sv.posedit.pos, hgvs.location.BaseOffsetInterval):
            return sv.posedit.pos.start.is_intronic or sv.posedit.pos.end.is_intronic
        return False

    def get_edit_type(self, sv: HgvsSequenceVariant) -> str | None:
        """Safely extract ``type`` property from SequenceVariant"""
        if sv is None or sv.posedit is None or sv.posedit.edit is None:
            return None
        return sv.posedit.edit.type

    def get_position_and_state(self, sv: HgvsSequenceVariant) -> tuple[int, int, str]:
        """Get the details of a sequence variant.

        Args:
            sv (hgvs.sequencevariant.SequenceVariant): The sequence variant object.

        Returns:
            tuple: A tuple containing the start position, end position, and state of the variant.

        Raises:
            ValueError: If the HGVS variant type is unsupported, or if ``pos``
                is an intronic :class:`~hgvs.location.BaseOffsetPosition`.
            TypeError: If ``pos`` is a :class:`~hgvs.location.BaseOffsetInterval`.

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
                state = self.data_proxy.get_sequence(sv.ac, start=start, end=end)
            else:
                state = sv.posedit.edit.alt or ""

        elif edit_type == "dup":
            ref = self.data_proxy.get_sequence(sv.ac, start=start, end=end)
            state = ref + ref

        else:
            msg = f"HGVS variant type {edit_type} is unsupported"
            raise ValueError(msg)

        return start, end, state

    def extract_allele_values(self, hgvs_expr: str) -> dict | None:
        """Parse hgvs into a VRS Allele Object

        kwargs:
            rle_seq_limit Optional(int): If RLE is set as the new state after
                normalization, this sets the limit for the length of the `sequence`.
                To exclude `sequence` from the response, set to 0.
                For no limit, set to `None`.
                Defaults value set in instance variable, `rle_seq_limit`.
            do_normalize (bool): `True` if fully justified normalization should be
                performed. `False` otherwise. Defaults to `True`

        #>>> a = tlr.from_hgvs("NC_000007.14:g.55181320A>T")
        #>>> a.model_dump()
        {
          'location': {
            'end': 55181320,
            'start': 55181319,
            'sequenceReference': {
              'type': 'SequenceReference',
              'refgetAccession': 'SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul'
            },
            'type': 'SequenceLocation'
          },
          'state': {
            'sequence': 'T',
            'type': 'LiteralSequenceExpression'
          },
          'type': 'Allele'
        }

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

    def from_allele(self, vo: models.Allele, namespace: str | None = None) -> list[str]:
        """Generate a *list* of HGVS expressions for VRS Allele.

        If `namespace` is not None, returns HGVS strings for the
        specified namespace.

        If `namespace` is None, returns HGVS strings for all alias
        translations.

        If no alias translations are available, an empty list is
        returned.

        If the VRS object cannot be expressed as HGVS, raises ValueError.

        This method assumes that IRIs are dereferenced, providing a `SequenceReference`
        as the `vo.location.sequenceReference`. If a `SequenceReference` is not
        provided, raises TypeError
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
        self, vo: models.Allele, sequence_type: str, sequence: str, accession: str
    ) -> HgvsSequenceVariant:
        """Create a SequenceVariant object from an Allele object."""
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
        """Perform hgvs library `normalize()` method"""
        return self.normalizer.normalize(hgvs)

    def n_to_c(self, hgvs: HgvsSequenceVariant) -> HgvsSequenceVariant:
        """Perform hgvs library `n_to_c` method"""
        return self.variant_mapper.n_to_c(hgvs)

    def c_to_n(self, hgvs: HgvsSequenceVariant) -> HgvsSequenceVariant:
        """Perform hgvs library `c_to_n` method"""
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
    hgvs_tools = HgvsTools(dp)
    hgvs_expr = hgvs_tools.from_allele(vrs_allele, namespace="refseq")

    print(hgvs_expr)  # noqa: T201
