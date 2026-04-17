"""Unit tests for the pure (non-networked) helpers in :mod:`ga4gh.vrs.utils.hgvs_tools`.

Uncertain-range HGVS expressions (e.g. ``(A_B)_(C_D)dup``, ``(?_N)_(M_?)del``)
require the hgvs library to be at ``2.0.0a0`` or newer, which exposes them as
nested :class:`hgvs.location.Interval` objects on ``sv.posedit.pos.start`` and
``.end``. See https://github.com/biocommons/hgvs/issues/225 for background.
"""

import hgvs.edit
import hgvs.exceptions
import hgvs.location
import hgvs.parser
import hgvs.posedit
import pytest

from ga4gh.vrs import models
from ga4gh.vrs.utils.hgvs_tools import _hgvs_pos_to_vrs, _vrs_pos_to_hgvs


@pytest.fixture(scope="module")
def hgvs_parser():
    return hgvs.parser.Parser()


class TestHgvsPosToVrs:
    """Tests for :func:`_hgvs_pos_to_vrs`."""

    def test_simple_position_start(self):
        pos = hgvs.location.SimplePosition(base=100)
        assert _hgvs_pos_to_vrs(pos, side="start") == 99

    def test_simple_position_end(self):
        pos = hgvs.location.SimplePosition(base=100)
        assert _hgvs_pos_to_vrs(pos, side="end") == 100

    def test_uncertain_range_start(self):
        pos = hgvs.location.Interval(
            start=hgvs.location.SimplePosition(base=100),
            end=hgvs.location.SimplePosition(base=200),
            uncertain=True,
        )
        result = _hgvs_pos_to_vrs(pos, side="start")
        assert isinstance(result, models.Range)
        assert result.root == [99, 199]

    def test_uncertain_range_end(self):
        pos = hgvs.location.Interval(
            start=hgvs.location.SimplePosition(base=100),
            end=hgvs.location.SimplePosition(base=200),
            uncertain=True,
        )
        result = _hgvs_pos_to_vrs(pos, side="end")
        assert isinstance(result, models.Range)
        assert result.root == [100, 200]

    def test_uncertain_range_with_unknown_lower_bound(self):
        pos = hgvs.location.Interval(
            start=hgvs.location.SimplePosition(base=None),
            end=hgvs.location.SimplePosition(base=200),
            uncertain=True,
        )
        result = _hgvs_pos_to_vrs(pos, side="start")
        assert isinstance(result, models.Range)
        assert result.root == [None, 199]

    def test_uncertain_range_with_unknown_upper_bound(self):
        pos = hgvs.location.Interval(
            start=hgvs.location.SimplePosition(base=100),
            end=hgvs.location.SimplePosition(base=None),
            uncertain=True,
        )
        result = _hgvs_pos_to_vrs(pos, side="end")
        assert isinstance(result, models.Range)
        assert result.root == [100, None]

    @pytest.mark.parametrize(
        ("hgvs_expr", "expected_start", "expected_end"),
        [
            # Balanced uncertain range (ClinVar 251062, GRCh37)
            (
                "NC_000019.9:g.(11211022_11213339)_(11217364_11218067)dup",
                [11211021, 11213338],
                [11217364, 11218067],
            ),
            # Left-unknown (ClinVar 425669)
            (
                "NC_000002.12:g.(?_202376935)_(202377551_202464808)del",
                [None, 202376934],
                [202377551, 202464808],
            ),
            # Right-unknown (ClinVar 425698)
            (
                "NC_000002.12:g.(202377551_202464808)_(202559947_?)del",
                [202377550, 202464807],
                [202559947, None],
            ),
            # Both sides unknown (ClinVar 220591)
            (
                "NC_000017.11:g.(?_58709859)_(58734342_?)del",
                [None, 58709858],
                [58734342, None],
            ),
        ],
    )
    def test_full_parse_to_vrs(
        self, hgvs_parser, hgvs_expr, expected_start, expected_end
    ):
        sv = hgvs_parser.parse_hgvs_variant(hgvs_expr)
        start = _hgvs_pos_to_vrs(sv.posedit.pos.start, side="start")
        end = _hgvs_pos_to_vrs(sv.posedit.pos.end, side="end")
        assert start.root == expected_start
        assert end.root == expected_end

    @pytest.mark.parametrize(
        ("hgvs_expr", "expected_start", "expected_end"),
        [
            # Exact left, uncertain right: start must be a plain int, not a
            # degenerate Range([99, 99]). Digest-level correctness depends on
            # this (a Range produces a different SequenceLocation digest).
            (
                "NC_000019.9:g.100_(200_?)del",
                99,
                models.Range(root=[200, None]),
            ),
            # Uncertain left, exact right: end must be a plain int.
            (
                "NC_000019.9:g.(?_200)_300del",
                models.Range(root=[None, 199]),
                300,
            ),
        ],
    )
    def test_mixed_endpoint_parse(
        self, hgvs_parser, hgvs_expr, expected_start, expected_end
    ):
        sv = hgvs_parser.parse_hgvs_variant(hgvs_expr)
        start = _hgvs_pos_to_vrs(sv.posedit.pos.start, side="start")
        end = _hgvs_pos_to_vrs(sv.posedit.pos.end, side="end")
        assert start == expected_start
        assert end == expected_end


class TestVrsPosToHgvs:
    """Tests for :func:`_vrs_pos_to_hgvs`."""

    def test_int_start(self):
        pos = _vrs_pos_to_hgvs(99, side="start")
        assert isinstance(pos, hgvs.location.SimplePosition)
        assert pos.base == 100

    def test_int_end(self):
        pos = _vrs_pos_to_hgvs(100, side="end")
        assert isinstance(pos, hgvs.location.SimplePosition)
        assert pos.base == 100

    def test_range_start(self):
        pos = _vrs_pos_to_hgvs(models.Range(root=[99, 199]), side="start")
        assert isinstance(pos, hgvs.location.Interval)
        assert pos.uncertain is True
        assert pos.start.base == 100
        assert pos.end.base == 200
        assert str(pos) == "(100_200)"

    def test_range_end(self):
        pos = _vrs_pos_to_hgvs(models.Range(root=[100, 200]), side="end")
        assert str(pos) == "(100_200)"

    def test_range_with_unknown_lower_bound(self):
        pos = _vrs_pos_to_hgvs(models.Range(root=[None, 199]), side="start")
        assert isinstance(pos, hgvs.location.Interval)
        assert str(pos) == "(?_200)"

    def test_range_with_unknown_upper_bound(self):
        pos = _vrs_pos_to_hgvs(models.Range(root=[100, None]), side="end")
        assert str(pos) == "(100_?)"

    def test_int_as_interval_wraps_in_non_uncertain_interval(self):
        pos = _vrs_pos_to_hgvs(99, side="start", as_interval=True)
        assert isinstance(pos, hgvs.location.Interval)
        assert pos.uncertain is False
        assert pos.start.base == 100
        assert pos.end.base == 100

    @pytest.mark.parametrize(
        "hgvs_expr",
        [
            "NC_000019.9:g.(11211022_11213339)_(11217364_11218067)dup",
            "NC_000002.12:g.(?_202376935)_(202377551_202464808)del",
            "NC_000017.11:g.(?_58709859)_(58734342_?)del",
        ],
    )
    def test_round_trip_position_strings(self, hgvs_parser, hgvs_expr):
        """Parsing an uncertain-range HGVS and feeding the positions back
        through the VRS <-> HGVS helpers should reproduce the original
        position strings.
        """
        sv = hgvs_parser.parse_hgvs_variant(hgvs_expr)
        start_vrs = _hgvs_pos_to_vrs(sv.posedit.pos.start, side="start")
        end_vrs = _hgvs_pos_to_vrs(sv.posedit.pos.end, side="end")
        start_hgvs = _vrs_pos_to_hgvs(start_vrs, side="start")
        end_hgvs = _vrs_pos_to_hgvs(end_vrs, side="end")
        assert str(start_hgvs) == str(sv.posedit.pos.start)
        assert str(end_hgvs) == str(sv.posedit.pos.end)


class TestUncertainRangePosEdit:
    """End-to-end formatting tests for the uncertain-range path in
    :meth:`HgvsTools._to_sequence_variant`. Network-free: exercises
    :func:`_vrs_pos_to_hgvs` together with the hgvs ``Interval`` +
    ``NARefAlt`` + ``PosEdit`` construction used there.
    """

    @pytest.mark.parametrize(
        ("start", "end", "expected"),
        [
            # Exact start, uncertain-range end.
            (99, models.Range(root=[200, 300]), "100_(200_300)del"),
            # Uncertain-range start, exact end.
            (models.Range(root=[99, 199]), 300, "(100_200)_300del"),
            # Uncertain on both sides.
            (
                models.Range(root=[99, 199]),
                models.Range(root=[300, 400]),
                "(100_200)_(300_400)del",
            ),
        ],
    )
    def test_mixed_and_full_uncertain_formatting(self, start, end, expected):
        ival = hgvs.location.Interval(
            start=_vrs_pos_to_hgvs(start, side="start", as_interval=True),
            end=_vrs_pos_to_hgvs(end, side="end", as_interval=True),
        )
        posedit = hgvs.posedit.PosEdit(
            pos=ival, edit=hgvs.edit.NARefAlt(ref="", alt=None)
        )
        assert str(posedit) == expected

    def test_ref_none_fails_for_uncertain_range(self):
        """The ``ref=""`` workaround in :meth:`HgvsTools._to_sequence_variant`
        exists because hgvs refuses to format a ``NARefAlt`` when both ``ref``
        and ``alt`` are undefined. This regression test locks that behavior
        in: if a future hgvs release relaxes the constraint, this test starts
        failing and the workaround can be re-evaluated.
        """
        ival = hgvs.location.Interval(
            start=_vrs_pos_to_hgvs(models.Range(root=[99, 199]), side="start"),
            end=_vrs_pos_to_hgvs(models.Range(root=[299, 399]), side="end"),
        )
        posedit = hgvs.posedit.PosEdit(
            pos=ival, edit=hgvs.edit.NARefAlt(ref=None, alt=None)
        )
        with pytest.raises(hgvs.exceptions.HGVSError):
            str(posedit)
