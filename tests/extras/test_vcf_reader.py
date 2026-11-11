"""Ensure proper functionality of VcfReader"""

import json
from pathlib import Path

import pytest

from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.annotator.vcf import VcfAnnotator, VcfReader, VcfReaderError

TEST_DATA_DIR = Path("tests/extras/data")


@pytest.fixture
def vcf_reader(dataproxy: SeqRepoDataProxy):
    return VcfReader(dataproxy)


@pytest.mark.parametrize("compute_for_ref", [True, False])
def test_read_vcf_round_trip(
    dataproxy: SeqRepoDataProxy,
    vcf_reader: VcfReader,
    tmp_path: Path,
    compute_for_ref: bool,
):
    """Test that the VRS Alleles read from an annotated VCF are the ones the annotator
    created
    """
    output_vcf = tmp_path / "test_vcf_output.vcf.gz"
    output_ndjson = tmp_path / "test_vcf_output.ndjson"
    VcfAnnotator(dataproxy).annotate(
        TEST_DATA_DIR / "test_vcf_input.vcf",
        output_vcf,
        vrs_attributes=True,
        compute_for_ref=compute_for_ref,
        output_ndjson_path=output_ndjson,
    )
    with output_ndjson.open() as f:
        expected_alleles = [json.loads(line) for line in f]

    read_alleles = []
    for record, vrs_alleles in vcf_reader.read(output_vcf):
        expected_num_alleles = (
            len(record.alleles) if compute_for_ref else len(record.alts)
        )
        assert len(vrs_alleles) == expected_num_alleles
        read_alleles += [
            vrs_allele.model_dump(exclude_none=True)
            for vrs_allele in vrs_alleles
            if vrs_allele is not None
        ]
    assert read_alleles == expected_alleles


def test_read_vcf(vcf_reader: VcfReader):
    """Test the VRS Alleles read for each kind of annotated record"""
    records = list(vcf_reader.read(TEST_DATA_DIR / "test_vcf_expected_output.vcf"))
    assert len(records) == 8
    vrs_alleles_by_pos = {record.pos: vrs_alleles for record, vrs_alleles in records}

    # REF is a ReferenceLengthExpression, ALT is a LiteralSequenceExpression
    ref, alt = vrs_alleles_by_pos[82664]
    assert ref.id == "ga4gh:VA.esmv8nARRRdSysDFuIErJxRAdUSVsWNE"
    assert ref.location.id == "ga4gh:SL.sXl66xA1nei_abGw9YmGkiAS-C9l3nUH"
    assert (
        ref.location.sequenceReference.refgetAccession
        == "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    )
    assert ref.location.start == 82663
    assert ref.location.end == 82664
    assert ref.state.model_dump(exclude_none=True) == {
        "type": "ReferenceLengthExpression",
        "length": 1,
        "sequence": "C",
        "repeatSubunitLength": 1,
    }
    assert alt.id == "ga4gh:VA.b9AIUegs1SOYO7nKY8l01yllkyUns-u4"
    assert alt.state.model_dump(exclude_none=True) == {
        "type": "LiteralSequenceExpression",
        "sequence": "T",
    }

    # Deletion whose empty ALT sequence is annotated as a missing VRS_States value
    _, alt = vrs_alleles_by_pos[490414]
    assert alt.id == "ga4gh:VA.lok7a3lot_cvUyw626otpJi4yxk0X07v"
    assert alt.state.model_dump(exclude_none=True) == {
        "type": "ReferenceLengthExpression",
        "length": 0,
        "sequence": "",
        "repeatSubunitLength": 2,
    }

    # Star allele
    ref, star, alt = vrs_alleles_by_pos[54220024]
    assert ref.id == "ga4gh:VA.LlmfhAC3gQlVQUwXWYiYjrn5V_K8vBz1"
    assert star is None
    assert alt.id == "ga4gh:VA.I7J3i1B36BACEUINcTwEh7uMv3I-PXT1"

    # Record that failed annotation
    assert vrs_alleles_by_pos[54220999] == [None, None]

    # ALT that couldn't be translated
    ref, alt, invalid_alt = vrs_alleles_by_pos[54221654]
    assert ref.id == "ga4gh:VA.kea5G-J1teg0iHMbgUELy-4L9lbJkgoj"
    assert alt.id == "ga4gh:VA.Zzlc24htmBV1HZZzWYgPD2_GfMInkrZu"
    assert invalid_alt is None


def test_read_vcf_altsonly(vcf_reader: VcfReader):
    """Test that VRS Alleles of a VCF annotated without REFs match its ALTs"""
    vrs_alleles_by_pos = {
        record.pos: vrs_alleles
        for record, vrs_alleles in vcf_reader.read(
            TEST_DATA_DIR / "test_vcf_expected_altsonly_output.vcf"
        )
    }
    (alt,) = vrs_alleles_by_pos[82664]
    assert alt.id == "ga4gh:VA.b9AIUegs1SOYO7nKY8l01yllkyUns-u4"
    star, alt = vrs_alleles_by_pos[54220024]
    assert star is None
    assert alt.id == "ga4gh:VA.I7J3i1B36BACEUINcTwEh7uMv3I-PXT1"
    assert vrs_alleles_by_pos[54220999] == [None]


@pytest.mark.parametrize(
    "input_vcf",
    ["test_vcf_input.vcf", "test_vcf_expected_output_no_vrs_attrs.vcf"],
)
def test_read_vcf_not_annotated(vcf_reader: VcfReader, input_vcf: str):
    """Test that a VCF without VRS attribute annotations is rejected"""
    with pytest.raises(VcfReaderError, match="missing VRS annotation INFO fields"):
        next(vcf_reader.read(TEST_DATA_DIR / input_vcf))


@pytest.mark.parametrize(
    ("annotation", "tampered_annotation", "error_match"),
    [
        (
            "VRS_Starts=82663,82663",
            "VRS_Starts=82662,82663",
            r"VRS Allele reconstructed at chr19-82664 has ID ga4gh:VA\.\S+, but was annotated with ga4gh:VA\.esmv8nARRRdSysDFuIErJxRAdUSVsWNE",
        ),
        (
            "VRS_Ends=82664,82664",
            "VRS_Ends=.,82664",
            "VRS Allele ga4gh:VA.esmv8nARRRdSysDFuIErJxRAdUSVsWNE at chr19-82664 is missing its VRS_Starts or VRS_Ends value",
        ),
        (
            "VRS_States=C,T",
            "VRS_States=C,t",
            "Invalid VRS annotations for VRS Allele ga4gh:VA.b9AIUegs1SOYO7nKY8l01yllkyUns-u4 at chr19-82664",
        ),
        (
            "VRS_Lengths=1,.",
            "VRS_Lengths=1",
            "Expected 2 VRS_Lengths values at chr19-82664, but found 1",
        ),
    ],
)
def test_read_vcf_invalid_annotations(
    vcf_reader: VcfReader,
    tmp_path: Path,
    annotation: str,
    tampered_annotation: str,
    error_match: str,
):
    """Test that a record whose VRS Alleles can't be reconstructed is rejected"""
    vcf_text = (TEST_DATA_DIR / "test_vcf_expected_output.vcf").read_text()
    assert annotation in vcf_text
    input_vcf = tmp_path / "test_vcf_tampered.vcf"
    input_vcf.write_text(vcf_text.replace(annotation, tampered_annotation, 1))

    with pytest.raises(VcfReaderError, match=error_match):
        list(vcf_reader.read(input_vcf))


@pytest.mark.parametrize(
    ("sequence", "length", "repeat_subunit_length", "expected"),
    [
        ("T", None, None, {"type": "LiteralSequenceExpression", "sequence": "T"}),
        (None, None, None, {"type": "LiteralSequenceExpression", "sequence": ""}),
        (
            "CTTTCTTT",
            8,
            4,
            {
                "type": "ReferenceLengthExpression",
                "length": 8,
                "sequence": "CTTTCTTT",
                "repeatSubunitLength": 4,
            },
        ),
        (
            None,
            0,
            2,
            {
                "type": "ReferenceLengthExpression",
                "length": 0,
                "sequence": "",
                "repeatSubunitLength": 2,
            },
        ),
        # sequence left out of the annotations because it's longer than RLE_SEQ_LIMIT
        (
            None,
            60,
            1,
            {
                "type": "ReferenceLengthExpression",
                "length": 60,
                "repeatSubunitLength": 1,
            },
        ),
        (None, 15, None, {"type": "LengthExpression", "length": 15}),
    ],
)
def test_get_state(sequence, length, repeat_subunit_length, expected):
    """Test that each kind of annotated state is reconstructed"""
    state = VcfReader._get_state(sequence, length, repeat_subunit_length)
    assert state.model_dump(exclude_none=True) == expected
