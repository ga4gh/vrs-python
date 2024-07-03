import pytest

from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import DataProxyValidationError
from ga4gh.vrs.extras.translator import AlleleTranslator


@pytest.fixture(scope="module")
def tlr(rest_dataproxy):
    return AlleleTranslator(
        data_proxy=rest_dataproxy,
        default_assembly_name="GRCh38",
        identify=False,
    )


# https://www.ncbi.nlm.nih.gov/clinvar/variation/17848/?new_evidence=true
snv_inputs = {
    "hgvs": "NC_000019.10:g.44908822C>T",
    "beacon": "19 : 44908822 C > T",
    "spdi": "NC_000019.10:44908821:1:T",
    "gnomad": "19-44908822-C-T"
}

snv_output = {
    "location": {
        "end": 44908822,
        "start": 44908821,
        "sequenceReference": {
            "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "T",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/693259/?new_evidence=true
mito_inputs = {
    "hgvs": "NC_012920.1:m.10083A>G",
    "beacon": "MT : 10083 A > G",
    "spdi": "NC_012920.1:10082:A:G",
    "gnomad": "MT-10083-A-G"
}

mito_output = {
    "location": {
      "start": 10082,
      "end": 10083,
      "sequenceReference": {
        "refgetAccession": "SQ.k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct",
        "type": "SequenceReference"
      },
      "type": "SequenceLocation"
    },
    "state": {
      "sequence": "G",
      "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1373966/?new_evidence=true
deletion_inputs = {
    "hgvs": "NC_000013.11:g.20003097del",
    "spdi": ["NC_000013.11:20003096:C:", "NC_000013.11:20003096:1:"],
    "gnomad": "13-20003096-AC-A"
}

deletion_output = {
    "location": {
        "end": 20003097,
        "start": 20003096,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

gnomad_deletion_output = {
    "location": {
        "end": 20003097,
        "start": 20003095,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "A",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

deletion_output_normalized = {
    "location": {
        "end": 20003097,
        "start": 20003096,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "length": 0,
        "repeatSubunitLength": 1,
        "sequence": "",
        "type": "ReferenceLengthExpression"
    },
    "type": "Allele"
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1687427/?new_evidence=true
insertion_inputs = {
    "hgvs": "NC_000013.11:g.20003010_20003011insG",
    "spdi": ["NC_000013.11:20003010::G", "NC_000013.11:20003010:0:G"],
    "gnomad": "13-20003010-A-AG"
}

insertion_output = {
    "location": {
        "end": 20003010,
        "start": 20003010,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "G",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

gnomad_insertion_output = {
    "location": {
        "end": 20003010,
        "start": 20003009,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "AG",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1264314/?new_evidence=true
duplication_inputs = {
    "hgvs": "NC_000013.11:g.19993838_19993839dup",
    "spdi": "NC_000013.11:19993837:GT:GTGT",
    "gnomad": "13-19993838-GT-GTGT"
}

duplication_output = {
    "location": {
        "end": 19993839,
        "start": 19993837,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "GTGT",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

duplication_output_normalized = {
    "location": {
        "end": 19993839,
        "start": 19993837,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "length": 4,
        "repeatSubunitLength": 2,
        "sequence": "GTGT",
        "type": "ReferenceLengthExpression"
    },
    "type": "Allele"
}

def test_from_invalid(tlr):
    try:
        tlr.translate_from("BRAF amplication")
        assert "Expected exception to be thrown" is None
    except ValueError as e:
        assert e.args[0] == "Unable to parse data as beacon, gnomad, hgvs, spdi, vrs"

@pytest.mark.vcr
def test_from_beacon(tlr):
    do_normalize = False
    assert tlr._from_beacon(snv_inputs["beacon"], do_normalize=do_normalize).model_dump(exclude_none=True) == snv_output
    assert tlr._from_beacon(mito_inputs["beacon"], do_normalize=do_normalize).model_dump(exclude_none=True) == mito_output


@pytest.mark.vcr
def test_from_gnomad(tlr):
    do_normalize = False
    assert tlr._from_gnomad(snv_inputs["gnomad"], do_normalize=do_normalize).model_dump(exclude_none=True) == snv_output
    assert tlr._from_gnomad(mito_inputs["gnomad"], do_normalize=do_normalize).model_dump(exclude_none=True) == mito_output
    assert tlr._from_gnomad(deletion_inputs["gnomad"], do_normalize=do_normalize).model_dump(exclude_none=True) == gnomad_deletion_output
    assert tlr._from_gnomad(insertion_inputs["gnomad"], do_normalize=do_normalize).model_dump(exclude_none=True) == gnomad_insertion_output
    assert tlr._from_gnomad(duplication_inputs["gnomad"], do_normalize=do_normalize).model_dump(exclude_none=True) == duplication_output

    # do_normalize defaults to true
    assert tlr._from_gnomad(snv_inputs["gnomad"]).model_dump(exclude_none=True) == snv_output
    assert tlr._from_gnomad(mito_inputs["gnomad"]).model_dump(exclude_none=True) == mito_output
    assert tlr._from_gnomad(deletion_inputs["gnomad"]).model_dump(exclude_none=True) == deletion_output_normalized
    assert tlr._from_gnomad(insertion_inputs["gnomad"]).model_dump(exclude_none=True) == insertion_output
    assert tlr._from_gnomad(duplication_inputs["gnomad"]).model_dump(exclude_none=True) == duplication_output_normalized

    assert tlr._from_gnomad("17-83129587-GTTGWCACATGA-G")

    # Test valid characters
    assert tlr._from_gnomad(
        "7-2-ACGTURYKMSWBDHVN-ACGTURYKMSWBDHVN",
        require_validation=False
    )

    # Invalid input. Ref does not match regex
    assert not tlr._from_gnomad("13-32936732-helloworld-C")

    # Ref != Actual ref
    invalid_var = "13-32936732-G-C"
    error_msg = "Expected reference sequence G on GRCh38:13 at positions (32936731, 32936732) but found C"

    with pytest.raises(DataProxyValidationError) as e:
        tlr._from_gnomad(invalid_var)
    assert str(e.value) == error_msg

    with pytest.raises(DataProxyValidationError) as e:
        tlr.translate_from(invalid_var, fmt="gnomad")
    assert str(e.value) == error_msg

    # require_validation set to False
    assert tlr._from_gnomad(invalid_var, require_validation=False)


@pytest.mark.vcr
def test_from_hgvs(tlr):
    do_normalize = False
    assert tlr._from_hgvs(snv_inputs["hgvs"], do_normalize=do_normalize).model_dump(exclude_none=True) == snv_output
    assert tlr._from_hgvs(mito_inputs["hgvs"], do_normalize=do_normalize).model_dump(exclude_none=True) == mito_output
    assert tlr._from_hgvs(deletion_inputs["hgvs"], do_normalize=do_normalize).model_dump(exclude_none=True) == deletion_output
    assert tlr._from_hgvs(insertion_inputs["hgvs"], do_normalize=do_normalize).model_dump(exclude_none=True) == insertion_output
    assert tlr._from_hgvs(duplication_inputs["hgvs"], do_normalize=do_normalize).model_dump(exclude_none=True) == duplication_output


@pytest.mark.vcr
def test_from_spdi(tlr):
    do_normalize = False
    assert tlr._from_spdi(snv_inputs["spdi"], do_normalize=do_normalize).model_dump(exclude_none=True) == snv_output
    assert tlr._from_spdi(mito_inputs["spdi"], do_normalize=do_normalize).model_dump(exclude_none=True) == mito_output
    for spdi_del_expr in deletion_inputs["spdi"]:
        assert tlr._from_spdi(spdi_del_expr, do_normalize=do_normalize).model_dump(exclude_none=True) == deletion_output, spdi_del_expr
    for spdi_ins_expr in insertion_inputs["spdi"]:
        assert tlr._from_spdi(spdi_ins_expr, do_normalize=do_normalize).model_dump(exclude_none=True) == insertion_output, spdi_ins_expr
    assert tlr._from_spdi(duplication_inputs["spdi"], do_normalize=do_normalize).model_dump(exclude_none=True) == duplication_output


@pytest.mark.vcr
def test_to_spdi(tlr):
    # do_normalize defaults to true
    spdiexpr = snv_inputs["spdi"]
    allele = tlr.translate_from(spdiexpr, "spdi")
    to_spdi = tlr.translate_to(allele, "spdi")
    assert 1 == len(to_spdi)
    assert spdiexpr == to_spdi[0]


hgvs_tests = (
    ("NC_000013.11:g.32936732=",
     {'digest': 'GJ2JySBMXePcV2yItyvCfbGBUoawOBON',
      'id': 'ga4gh:VA.GJ2JySBMXePcV2yItyvCfbGBUoawOBON',
      'location': {'digest': '28YsnRvD40gKu1x3nev0gRzRz-5OTlpS',
                   'end': 32936732,
                   'id': 'ga4gh:SL.28YsnRvD40gKu1x3nev0gRzRz-5OTlpS',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 32936731,
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'C', 'type': 'LiteralSequenceExpression'},
      'type': 'Allele'}),
    ("NC_000007.14:g.55181320A>T",
     {'digest': 'Hy2XU_-rp4IMh6I_1NXNecBo8Qx8n0oE',
      'id': 'ga4gh:VA.Hy2XU_-rp4IMh6I_1NXNecBo8Qx8n0oE',
      'location': {'digest': '_G2K0qSioM74l_u3OaKR0mgLYdeTL7Xd',
                   'end': 55181320,
                   'id': 'ga4gh:SL._G2K0qSioM74l_u3OaKR0mgLYdeTL7Xd',
                   'sequenceReference': {'refgetAccession': 'SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                                         'type': 'SequenceReference'},
                   'start': 55181319,
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'T', 'type': 'LiteralSequenceExpression'},
      'type': 'Allele'}),
    ("NC_000007.14:g.55181220del",
     {'digest': 'klRMVChjvV73ZxS9Ajq1Rb8WU-p_HbLu',
      'id': 'ga4gh:VA.klRMVChjvV73ZxS9Ajq1Rb8WU-p_HbLu',
      'location': {'digest': 'ljan7F0ePe9uiD6f2u80ZG5gDtx9Mr0V',
                   'end': 55181220,
                   'id': 'ga4gh:SL.ljan7F0ePe9uiD6f2u80ZG5gDtx9Mr0V',
                   'sequenceReference': {'refgetAccession': 'SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                                         'type': 'SequenceReference'},
                   'start': 55181219,
                   'type': 'SequenceLocation'},
      'state': {'length': 0,
                'repeatSubunitLength': 1,
                'sequence': '',
                'type': 'ReferenceLengthExpression'},
      'type': 'Allele'}),
    ("NC_000007.14:g.55181230_55181231insGGCT",
     {'digest': 'CLOvnFRJXGNRB9aTuNbvsLqc7syRYb55',
      'id': 'ga4gh:VA.CLOvnFRJXGNRB9aTuNbvsLqc7syRYb55',
      'location': {'digest': 'lh4dRt_xWPi3wrubcfomi5DkD7fu6wd2',
                   'end': 55181230,
                   'id': 'ga4gh:SL.lh4dRt_xWPi3wrubcfomi5DkD7fu6wd2',
                   'sequenceReference': {'refgetAccession': 'SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                                         'type': 'SequenceReference'},
                   'start': 55181230,
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'GGCT', 'type': 'LiteralSequenceExpression'},
      'type': 'Allele'}),
    ("NC_000013.11:g.32331093_32331094dup",
     {'digest': 'swY2caCgv1kP6YqKyPlcEzJqTvou15vC',
      'id': 'ga4gh:VA.swY2caCgv1kP6YqKyPlcEzJqTvou15vC',
      'location': {'digest': 'ikECYncPpE1xh6f_LiComrFGevocjDHQ',
                   'end': 32331094,
                   'id': 'ga4gh:SL.ikECYncPpE1xh6f_LiComrFGevocjDHQ',
                   'sequenceReference': {
                       'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                       'type': 'SequenceReference'},
                   'start': 32331082,
                   'type': 'SequenceLocation'},
      'state': {'length': 14,
                'repeatSubunitLength': 2,
                'sequence': 'TTTTTTTTTTTTTT',
                'type': 'ReferenceLengthExpression'},
      'type': 'Allele'}),
    ("NC_000013.11:g.32316467dup",
     {'digest': '96ak7XdY3DNbp71aHEXw-NHSfeHGW-KT',
      'id': 'ga4gh:VA.96ak7XdY3DNbp71aHEXw-NHSfeHGW-KT',
      'location': {'digest': 'fwfHu8VaD2-6Qvay9MJSINXPS767RYSw',
                   'end': 32316467,
                   'id': 'ga4gh:SL.fwfHu8VaD2-6Qvay9MJSINXPS767RYSw',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 32316466,
                   'type': 'SequenceLocation'},
      'state': {'length': 2,
                'repeatSubunitLength': 1,
                'sequence': 'AA',
                'type': 'ReferenceLengthExpression'},
      'type': 'Allele'}),
    ("NM_001331029.1:c.722A>G",
     {'digest': 'DPe4AO-S0Yu4wzSCmys7eGn4p4sO0zaC',
      'id': 'ga4gh:VA.DPe4AO-S0Yu4wzSCmys7eGn4p4sO0zaC',
      'location': {'digest': '7hcVmPnIspQNDfZKBzRJFc8K9GaJuAlY',
                   'end': 872,
                   'id': 'ga4gh:SL.7hcVmPnIspQNDfZKBzRJFc8K9GaJuAlY',
                   'sequenceReference': {'refgetAccession': 'SQ.MBIgVnoHFw34aFqNUVGM0zgjC3d-v8dK',
                                         'type': 'SequenceReference'},
                   'start': 871,
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'G', 'type': 'LiteralSequenceExpression'},
      'type': 'Allele'}),
    ("NM_181798.1:c.1007G>T",
     {'digest': 'vSL4aV7mPQKQLX7Jk-PmXN0APs0cBIr9',
      'id': 'ga4gh:VA.vSL4aV7mPQKQLX7Jk-PmXN0APs0cBIr9',
      'location': {'digest': 'EtvHvoj1Lsq-RruzIzWbKOIAW-bt193w',
                   'end': 1263,
                   'id': 'ga4gh:SL.EtvHvoj1Lsq-RruzIzWbKOIAW-bt193w',
                   'sequenceReference': {'refgetAccession': 'SQ.KN07u-RFqd1dTyOWOG98HnOq87Nq-ZIg',
                                         'type': 'SequenceReference'},
                   'start': 1262,
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'T', 'type': 'LiteralSequenceExpression'},
      'type': 'Allele'}),
    ("NC_000019.10:g.289464_289465insCACA",
     {'digest': 'YFUR4oR_84b-rRFf0UzOjfI4eE5FTKAP',
      'id': 'ga4gh:VA.YFUR4oR_84b-rRFf0UzOjfI4eE5FTKAP',
      'type': 'Allele',
      'location': {'digest': 'L145KFLJeJ334YnOVm59pPlbdqfHhgXZ',
                   'end': 289466,
                   'id': 'ga4gh:SL.L145KFLJeJ334YnOVm59pPlbdqfHhgXZ',
                   'sequenceReference': {'refgetAccession': 'SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl',
                                         'type': 'SequenceReference'},
                   'start': 289464,
                   'type': 'SequenceLocation'},
        'state': {'length': 6,
                  'repeatSubunitLength': 2,
                  'sequence': 'CACACA',
                  'type': 'ReferenceLengthExpression'}}),
    ("NC_000019.10:g.289485_289500del",
     {'digest': 'Djc_SwVDFunsArqwUM00PciVaF70VTcU',
      'id': 'ga4gh:VA.Djc_SwVDFunsArqwUM00PciVaF70VTcU',
      'type': 'Allele',
      'location': {'digest': 'WTE7jyihK4qvRRzEqM7u5nSD4iS2k3xp',
                   'end': 289501,
                   'id': 'ga4gh:SL.WTE7jyihK4qvRRzEqM7u5nSD4iS2k3xp',
                   'sequenceReference': {'refgetAccession': 'SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl',
                                         'type': 'SequenceReference'},
                   'start': 289480,
                   'type': 'SequenceLocation'},
        'state': {'length': 5,
                  'repeatSubunitLength': 16,
                  'sequence': 'CGAGG',
                  'type': 'ReferenceLengthExpression'}}),
)

hgvs_tests_to_hgvs_map = {
    "NC_000019.10:g.289464_289465insCACA": "NC_000019.10:g.289466_289467insCACA",
    "NC_000019.10:g.289485_289500del": "NC_000019.10:g.289486_289501del"
}


@pytest.mark.parametrize("hgvsexpr,expected", hgvs_tests)
@pytest.mark.vcr
def test_hgvs(tlr, hgvsexpr, expected):
    # do_normalize defaults to true
    tlr.identify = True
    allele = tlr.translate_from(hgvsexpr, "hgvs")
    assert allele.model_dump(exclude_none=True) == expected

    to_hgvs = tlr.translate_to(allele, "hgvs")
    assert (hgvsexpr in to_hgvs) or (hgvs_tests_to_hgvs_map.get(hgvsexpr, hgvsexpr) in to_hgvs)

@pytest.mark.vcr
def test_rle_seq_limit(tlr):
    # do_normalize defaults to true
    tlr.identify = True

    a_dict = {
        'digest': 'j7qUzb1uvmdxLAbtdCPiay4kIRQmyZNv',
        'id': 'ga4gh:VA.j7qUzb1uvmdxLAbtdCPiay4kIRQmyZNv',
        'location': {'digest': '88oOqkUgALP7fnN8P8lbvCosFhG8YpY0',
                     'end': 32331094,
                     'id': 'ga4gh:SL.88oOqkUgALP7fnN8P8lbvCosFhG8YpY0',
                     'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                           'type': 'SequenceReference'},
                     'start': 32331042,
                     'type': 'SequenceLocation'},
        'state': {'length': 104,
                  'repeatSubunitLength': 52,
                  'type': 'ReferenceLengthExpression'},
        'type': 'Allele'}
    input_hgvs_expr = "NC_000013.11:g.32331043_32331094dup"

    # use default rle_seq_limit
    allele_no_seq = tlr.translate_from(input_hgvs_expr, fmt="hgvs")
    assert allele_no_seq.model_dump(exclude_none=True) == a_dict

    with pytest.raises(AttributeError, match="'NoneType' object has no attribute 'root'"):
        tlr.translate_to(allele_no_seq, "hgvs")

    # set rle_seq_limit to None
    allele_with_seq = tlr.translate_from(input_hgvs_expr, fmt="hgvs", rle_seq_limit=None)
    a_dict_with_seq = a_dict.copy()
    a_dict_with_seq["state"][
        "sequence"] = "TTTAGTTGAACTACAGGTTTTTTTGTTGTTGTTGTTTTGATTTTTTTTTTTTTTTAGTTGAACTACAGGTTTTTTTGTTGTTGTTGTTTTGATTTTTTTTTTTT"
    assert allele_with_seq.model_dump(exclude_none=True) == a_dict_with_seq

    output_hgvs_expr = tlr.translate_to(allele_with_seq, "hgvs")
    assert output_hgvs_expr == [input_hgvs_expr]

@pytest.mark.vcr
def test_to_hgvs_iri_ref_keyerror(tlr):
    # IRI is passed
    iri_vo = models.Allele(
        **{
            "location": {
                "end": 1263,
                "start": 1262,
                "sequenceReference": "seqrefs.jsonc#/NM_181798.1",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        }
    )
    with pytest.raises(KeyError) as e:
        # even though the seqrefs.jsonc#/NM_181798.1 is a valid iri-reference for json schema, it is not yet handled here
        # we have to add functionality to address the handling of iri-references in the future
        tlr.translate_to(iri_vo, "hgvs")
    assert str(e.value) == "'ga4gh:seqrefs.jsonc#/NM_181798.1'"

# TODO: Readd these tests
# @pytest.mark.vcr
# def test_errors(tlr):
#     with pytest.raises(ValueError):
#         tlr._from_beacon("bogus")
#
#     with pytest.raises(ValueError):
#         tlr._from_gnomad("NM_182763.2:c.688+403C>T")
#
#     with pytest.raises(ValueError):
#         tlr._from_hgvs("NM_182763.2:c.688+403C>T")
#
#     with pytest.raises(ValueError):
#         tlr._from_hgvs("NM_182763.2:c.688_690inv")
#
#     with pytest.raises(ValueError):
#         tlr._from_spdi("NM_182763.2:c.688+403C>T")
