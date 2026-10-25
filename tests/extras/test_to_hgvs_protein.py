"""Tests for to_hgvs support of protein variants (single-residue changes).

See https://github.com/ga4gh/vrs-python/issues/633.

These tests use a stub data proxy and a mocked UTA connection so they run
without SeqRepo data or external services.
"""

from unittest import mock

import hgvs.exceptions
import pytest

from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import _DataProxy
from ga4gh.vrs.utils.hgvs_tools import HgvsTools


class _StubProteinDataProxy(_DataProxy):
    """Minimal data proxy stub exposing a single protein sequence.

    Residue 261 (1-based) of refseq:NP_060204.1 is Val.
    """

    REFGET_ACCESSION = "SQ.dvlWjX2CGulwfb2ehmkCFn02ah7tEEVB"
    ACCESSION = "NP_060204.1"

    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        assert identifier == f"ga4gh:{self.REFGET_ACCESSION}"
        assert (start, end) == (260, 261)
        return "V"

    def get_metadata(self, identifier: str) -> dict:
        assert identifier == f"ga4gh:{self.REFGET_ACCESSION}"
        return {
            "aliases": [
                f"ga4gh:{self.REFGET_ACCESSION}",
                f"refseq:{self.ACCESSION}",
            ],
            "length": 597,
            "alphabet": "ACDEFGHIKLMNPQRSTVWXY",
        }


@pytest.fixture(scope="module")
def hgvs_tools():
    """HgvsTools with stubbed data proxy and no UTA connection.

    ``normalize`` is stubbed to raise HGVSDataNotAvailableError so the
    no-UTA-data fallback path is exercised deterministically.
    """
    with (
        mock.patch("hgvs.dataproviders.uta.connect", return_value=None),
        mock.patch.object(
            HgvsTools,
            "normalize",
            side_effect=hgvs.exceptions.HGVSDataNotAvailableError("no UTA data"),
        ),
    ):
        yield HgvsTools(_StubProteinDataProxy())


def _protein_allele(
    start: int = 260,
    end: int = 261,
    state: dict | None = None,
) -> models.Allele:
    return models.Allele.model_validate(
        {
            "type": "Allele",
            "location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": _StubProteinDataProxy.REFGET_ACCESSION,
                },
                "start": start,
                "end": end,
            },
            "state": state or {"type": "LiteralSequenceExpression", "sequence": "A"},
        }
    )


def test_to_hgvs_protein_substitution(hgvs_tools):
    """The exact example from issue #633."""
    allele = models.Allele.model_validate(
        {
            "id": "ga4gh:VA.AIm-GH_iqp_bpcIVzi431fl7z-cQimNN",
            "type": "Allele",
            "digest": "AIm-GH_iqp_bpcIVzi431fl7z-cQimNN",
            "location": {
                "id": "ga4gh:SL.wgFs8Z2Nk4uD7nq_gFbdYH1FhagXkMhL",
                "type": "SequenceLocation",
                "digest": "wgFs8Z2Nk4uD7nq_gFbdYH1FhagXkMhL",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": _StubProteinDataProxy.REFGET_ACCESSION,
                },
                "start": 260,
                "end": 261,
            },
            "state": {"type": "LiteralSequenceExpression", "sequence": "A"},
        }
    )
    assert hgvs_tools.from_allele(allele, namespace="refseq") == [
        "NP_060204.1:p.Val261Ala"
    ]


def test_to_hgvs_protein_nonsense(hgvs_tools):
    """Single-residue nonsense substitution -> p.Val261Ter."""
    allele = _protein_allele(
        state={"type": "LiteralSequenceExpression", "sequence": "*"}
    )
    assert hgvs_tools.from_allele(allele, namespace="refseq") == [
        "NP_060204.1:p.Val261Ter"
    ]


def test_to_hgvs_protein_deletion(hgvs_tools):
    """Single-residue deletion (empty alt) -> p.Val261del."""
    allele = _protein_allele(
        state={"type": "LiteralSequenceExpression", "sequence": ""}
    )
    assert hgvs_tools.from_allele(allele, namespace="refseq") == [
        "NP_060204.1:p.Val261del"
    ]


def test_to_hgvs_protein_multi_residue_rejected(hgvs_tools):
    """Alleles spanning more than one residue are not supported (yet)."""
    allele = _protein_allele(start=260, end=262)
    with pytest.raises(ValueError, match="single-residue"):
        hgvs_tools.from_allele(allele, namespace="refseq")


def test_to_hgvs_protein_insertion_rejected(hgvs_tools):
    """Insertions (start == end) are not single-residue changes."""
    allele = _protein_allele(start=260, end=260)
    with pytest.raises(ValueError, match="single-residue"):
        hgvs_tools.from_allele(allele, namespace="refseq")


def test_to_hgvs_protein_rle_state_rejected(hgvs_tools):
    """ReferenceLengthExpression states are not supported for proteins."""
    allele = _protein_allele(
        state={
            "type": "ReferenceLengthExpression",
            "length": 1,
            "repeatSubunitLength": 1,
            "sequence": "V",
        }
    )
    with pytest.raises(ValueError, match="LiteralSequenceExpression"):
        hgvs_tools.from_allele(allele, namespace="refseq")


def test_to_hgvs_protein_reference_allele_rejected(hgvs_tools):
    """An allele identical to reference cannot be expressed as a variant."""
    allele = _protein_allele(
        state={"type": "LiteralSequenceExpression", "sequence": "V"}
    )
    with pytest.raises(ValueError, match="Reference alleles"):
        hgvs_tools.from_allele(allele, namespace="refseq")


@pytest.mark.parametrize("alt", ["AA"])
def test_to_hgvs_protein_multi_residue_alt_rejected(hgvs_tools, alt):
    """Alternate alleles longer than one residue are rejected."""
    allele = _protein_allele(
        state={"type": "LiteralSequenceExpression", "sequence": alt}
    )
    with pytest.raises(ValueError, match="Unsupported protein"):
        hgvs_tools.from_allele(allele, namespace="refseq")


def test_to_hgvs_protein_invalid_alt_rejected(hgvs_tools):
    """Alternate residues that are not valid amino acid codes are rejected.

    (Constructed via attribute assignment since the VRS model itself
    rejects such sequences at validation time.)
    """
    allele = _protein_allele()
    allele.state.sequence.root = "!"
    with pytest.raises(ValueError, match="Unsupported protein"):
        hgvs_tools.from_allele(allele, namespace="refseq")
