import pytest
from tempfile import TemporaryFile
from ga4gh.vrs import models

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


# provide test inputs and expected outputs testing Translator._from_vcf_record
to_vcf_tests = (
    # SNPs
    (
        ('1', '92633', 'C', ['T']),
        {
            '_id': 'ga4gh:VA.qAK6JCN3-AVa9_6Qq3AqAuppUU0bWgfH',
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 92632,
                    'end': 92633
                }
            },
            'state': {'type': 'SequenceState', 'sequence': 'T'}
        }
    ),
    (
        ('1', '63002', 'A', ['G']),
        {
           '_id': 'ga4gh:VA.VewFnlxS7DmjEdKMkj0xZK-9GaHGMJcx',
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO',
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
            '_id': 'ga4gh:VA.TDba99qQKLMUiooxOgIX_EEFDgyNPBAq',
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 22304601,
                    'end': 22304607
                }
            },
            'state': {'type': 'SequenceState', 'sequence': 'AAAAAAA'}
        }
    ),
    (
        ('1', '72297', 'G', ['GTAT']),
        {
            '_id': 'ga4gh:VA.wbhpDCQ0MRtG0pZZWH-yarqSdOjGNJEL',
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO',
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
            '_id': 'ga4gh:VA.mId9FgDqwHkP3mBt1lgfSr2-bhCJaxGg',
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 29204173, 'end': 29204179
                }
            },
            'state': {'type': 'SequenceState', 'sequence': 'ACA'}
        }
    ),
    (
        ('4', '116619313', 'GT', ['G']),
        {
            '_id': 'ga4gh:VA.GVcxxtyjvhuuCtcdpcNAvVhRHWbzAhra',
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:SQ.HxuclGHh0XCDuF8x6yQrpHUBL7ZntAHc',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 116619313,
                    'end': 116619314
                }
            },
            'state': {'type': 'SequenceState', 'sequence': ''}
        }
    )
)


@pytest.mark.parametrize("record,expected", to_vcf_tests)
def test_from_vcf_record(tlr_norm, record, expected):
    """Test Translator._from_vcf_record"""
    tlr_norm.normalize = True
    allele = tlr_norm._from_vcf_record(*record)
    assert allele.as_dict() == expected


@pytest.fixture(scope="session")
def vcf_file_in():
    """Provide input file fixture to test Translator._from_vcf"""
    file_string = '##fileformat=VCFv4.3\n' + \
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
                  '20\t14370\trs6054257\tG\tA\t29\tPASS\tDP=14;AF=0.5;DB\tGT:DP\t0/0:1\t0/1:8\t1/1:5\n' + \
                  '4\t116612759\trs566366173\tATTGTT\tA\t.\tPASS\tRS=566366173;RSPOS=116612760;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000000005040026000200;WGT=1;VC=DIV;ASP;VLD;KGPhase3;CAF=0.9968,0.003195;COMMON=1;TOPMED=0.99684633027522935,0.00315366972477064'
    file_bytes = file_string.encode('utf-8')
    fp = TemporaryFile()
    fp.write(file_bytes)
    fp.seek(0)
    return fp


@pytest.fixture(scope="session")
def vcf_from_file_expected():
    """Provide expected output to test Translator._from_vcf"""
    return [
        {
            '_id': 'ga4gh:VA.m8e1aRLy--eKkjzCJWz5w7fZEBmOhbXZ',
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:SQ.-A1QmD_MatoqxvgVxBLZTONHz9-c7nQo',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 14369,
                    'end': 14370
                }
            },
            'state': {'type': 'SequenceState', 'sequence': 'A'}
        },
        {
            '_id': 'ga4gh:VA.EG07sjR0AIU2XOTHvlXvM__egaimtuPT',
            'type': 'Allele',
            'location': {
                'type': 'SequenceLocation',
                'sequence_id': 'ga4gh:SQ.HxuclGHh0XCDuF8x6yQrpHUBL7ZntAHc',
                'interval': {
                    'type': 'SimpleInterval',
                    'start': 116612759,
                    'end': 116612766
                }
            },
            'state': {'type': 'SequenceState', 'sequence': 'TT'}
        }
    ]


def test_from_vcf_file(tlr_norm, vcf_file_in, vcf_from_file_expected):
    """Test Translator._from_vcf"""
    alleles = tlr_norm._from_vcf(vcf_file_in)
    for i in range(len(alleles)):
        assert alleles[i].as_dict() == vcf_from_file_expected[i]


@pytest.fixture(scope="session")
def vcf_file_vrs_input():
    """Provide VRS Allele objects as input to test Translator._to_vcf"""
    return [
        models.Allele(
            _id='ga4gh:VA.VewFnlxS7DmjEdKMkj0xZK-9GaHGMJcx',
            location=models.Location(
                sequence_id='ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO',
                interval=models.SimpleInterval(start=63001, end=63002)
            ),
            state=models.SequenceState(sequence='G')
        ),
        models.Allele(
            _id='ga4gh:VA.qAK6JCN3-AVa9_6Qq3AqAuppUU0bWgfH',
            location=models.Location(
                sequence_id='ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO',
                interval=models.SimpleInterval(start=92632, end=92633)
            ),
            state=models.SequenceState(sequence='T')
        )
    ]

@pytest.fixture(scope="session")
def vcf_from_file_expected():
    """Provide expected output for Translator._to_vcf"""
    return [
        '1\t63002\t.\tA\tG\t.\t.\t.\n',
        '1\t92633\t.\tC\tT\t.\t.\t.\n'
    ]

def test_to_vcf(tlr_norm, vcf_file_vrs_input, vcf_from_file_expected):
    """Test Translator._to_vcf"""
    outfile_path = "test_out.vcf"
    tlr_norm._to_vcf(vcf_file_vrs_input, outfile_path)
    with open(outfile_path) as f:
        outfile_lines = list(f.readlines())
    assert outfile_lines[4] == vcf_from_file_expected[0]
    assert outfile_lines[5] == vcf_from_file_expected[1]
