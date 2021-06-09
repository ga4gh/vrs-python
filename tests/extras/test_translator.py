import pytest
from tempfile import TemporaryFile

inputs = {
    "hgvs": "NC_000013.11:g.32936732G>C",
    "beacon": "13 : 32936732 G > C",
    "spdi": "NC_000013.11:32936731:1:C",
    "gnomad": "13-32936732-G-C",
    "vcf": ['13', '32936732', 'G', 'C']
}

output = {
    'location': {
        'interval': {'end': 32936732, 'start': 32936731, 'type': 'SimpleInterval'},
        'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
        'type': 'SequenceLocation'},
    'state': {'sequence': 'C', 'type': 'SequenceState'},
    'type': 'Allele'
}


@pytest.mark.vcr
def test_from_beacon(tlr):
    assert tlr._from_beacon(inputs["beacon"]).as_dict() == output


@pytest.mark.vcr
def test_from_gnomad(tlr):
    assert tlr._from_gnomad(inputs["gnomad"]).as_dict() == output


@pytest.mark.vcr
def test_from_hgvs(tlr):
    assert tlr._from_hgvs(inputs["hgvs"]).as_dict() == output


@pytest.mark.vcr
def test_from_spdi(tlr):
    assert tlr._from_spdi(inputs["spdi"]).as_dict() == output


@pytest.mark.vcr
def test_from_vcf(tlr):
    assert tlr._from_vcf(inputs["vcf"]).as_dict() == output

hgvs_tests = (
    ("NC_000013.11:g.32936732=",
     {'location': {'interval': {'end': 32936732,
                                'start': 32936731,
                                'type': 'SimpleInterval'},
                   'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'C', 'type': 'SequenceState'},
      'type': 'Allele'}),

    ("NC_000007.14:g.55181320A>T",
     {'location': {'interval': {'end': 55181320,
                                'start': 55181319,
                                'type': 'SimpleInterval'},
                   'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'T', 'type': 'SequenceState'},
      'type': 'Allele'}),

    ("NC_000007.14:g.55181220del",
     {'location': {'interval': {'end': 55181220,
                                'start': 55181219,
                                'type': 'SimpleInterval'},
                   'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                   'type': 'SequenceLocation'},
      'state': {'sequence': '', 'type': 'SequenceState'},
      'type': 'Allele'}),

    ("NC_000007.14:g.55181230_55181231insGGCT",
     {'location': {'interval': {'end': 55181230,
                                'start': 55181230,
                                'type': 'SimpleInterval'},
                   'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'GGCT', 'type': 'SequenceState'},
      'type': 'Allele'}),
    ("NC_000013.11:g.32331093_32331094dup",
     {'location': {'interval': {'end': 32331094,
                                'start': 32331082, 'type': 'SimpleInterval'},
                   'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'TTTTTTTTTTTTTT', 'type': 'SequenceState'},
      'type': 'Allele'}),
     ("NC_000013.11:g.32316467dup",
      {'location': {'interval': {'end': 32316467,
                                'start': 32316466,
                                'type': 'SimpleInterval'},
                   'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'AA', 'type': 'SequenceState'},
       'type': 'Allele'}),
)


@pytest.mark.parametrize("hgvsexpr,expected", hgvs_tests)
@pytest.mark.vcr
def test_hgvs(tlr, hgvsexpr, expected):
    tlr.normalize = True
    allele = tlr.translate_from(hgvsexpr, "hgvs")
    assert expected == allele.as_dict()

    to_hgvs = tlr.translate_to(allele, "hgvs")
    assert 1 == len(to_hgvs)
    assert hgvsexpr == to_hgvs[0]

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


vcf_tests = (
    # SNPs
    (
        ('1', '92633', 'C', ['T']),
        {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:GS.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 92632, 'end': 92633
                }
            },
            'state': {'type': 'SequenceState', 'sequence': 'T'}
        }
    ),
    (
        ('1', '63002', 'A', ['G']),
        {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:GS.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 63001,
                    'end': 63002
                }
            },
            'state': {'type': 'SequenceState', 'sequence': 'G'}
        }
    ),
    # ins
    (
        ('Y', '22304601', 'G', ['GA']),
        {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:GS.BT7QyW5iXaX_1PSX-msSGYsqRdMKqkj-',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 22304601,
                    'end': 22304601
                }
            },
            'state': {'type': 'SequenceState', 'sequence': 'A'}
        }
    ),
    (
        ('1', '72297', 'G', ['GTAT']),
        {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:GS.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 72297,
                    'end': 72301
                }
            },
            'state': {'type': 'SequenceState', 'sequence': 'TATTATT'}
        }
    ),
    # del
    (
        ('17', '29204173', 'TACA', ['T']),
        {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:GS.AjWXsI7AkTK35XW9pgd3UbjpC3MAevlz',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 29204173,
                    'end': 29204176
                }
            },
            'state': {'type': 'SequenceState', 'sequence': ''}
        }
    ),
    (
        ('4', '116619313', 'GT', ['G']),
        {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:GS.iy7Zfceb5_VGtTQzJ-v5JpPbpeifHD_V',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 116619312,
                    'end': 116619314
                }
            },
            'state': {'type': 'SequenceState', 'sequence': 'G'}
        }
    )
)


@pytest.mark.parametrize("record,expected", vcf_tests)
def test_vcf(tlr, record, expected):
    tlr.normalize = True
    allele = tlr._from_vcf_record(*record, 'GRCh37')
    assert allele.as_dict() == expected


vcf_files = (
    [
        '##fileformat=VCFv4.3\n' + \
        '##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta\n' + \
        '##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>\n' + \
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n' + \
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n' + \
        '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">\n' + \
        '##FILTER=<ID=q10,Description="Quality below 10">\n' + \
        '##FILTER=<ID=s50,Description="Less than 50% of samples have data">\n' + \
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' + \
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n' + \
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003\n' + \
        '20\t14370\trs6054257\tG\tA\t29\tPASS\tDP=14;AF=0.5;DB\tGT:DP\t0/0:1\t0/1:8\t1/1:5\n',
        {
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:GS.-A1QmD_MatoqxvgVxBLZTONHz9-c7nQo',
                'interval': {
                    'type': 'SimpleInterval', 'start': 14369, 'end': 14370
                }
            },
            'state': {
                'type': 'SequenceState', 'sequence': 'A'
            }
        }
    ],
)




@pytest.mark.parametrize("vcf_file,vcf_file_expected", vcf_files)
def test_vcf_file(tlr, vcf_file, vcf_file_expected):
    v = vcf_file.encode('utf-8')
    fp = TemporaryFile()
    fp.write(v)
    fp.seek(0)
    alleles = tlr._from_vcf(fp)
    assert alleles[0].as_dict() == vcf_file_expected
