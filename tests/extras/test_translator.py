import pytest
from tempfile import TemporaryFile
from ga4gh.vrs import models
from os import remove
import time
from pysam import VariantHeader

inputs = {
    "hgvs": "NC_000013.11:g.32936732G>C",
    "beacon": "13 : 32936732 G > C",
    "spdi": "NC_000013.11:32936731:1:C",
    "gnomad": "13-32936732-G-C",
    "vcf": ['13', '32936732', 'G', ['C']]
}

output = {
    'location': {
        'interval': {
            'end': 32936732,
            'start': 32936731,
            'type': 'SimpleInterval'
        },
        'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
        'type': 'SequenceLocation'
    },
    'state': {
        'sequence': 'C',
        'type': 'SequenceState'
    },
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
    assert tlr._from_vcf_record(*inputs["vcf"], assembly_name='GRCh38')[0].as_dict() == output


hgvs_tests = (
    ("NC_000013.11:g.32936732=", {
        'location': {
            'interval': {
                'end': 32936732,
                'start': 32936731,
                'type': 'SimpleInterval'
            },
            'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'C',
            'type': 'SequenceState'
        },
        'type': 'Allele'
    }),
    ("NC_000007.14:g.55181320A>T", {
        'location': {
            'interval': {
                'end': 55181320,
                'start': 55181319,
                'type': 'SimpleInterval'
            },
            'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'T',
            'type': 'SequenceState'
        },
        'type': 'Allele'
    }),
    ("NC_000007.14:g.55181220del", {
        'location': {
            'interval': {
                'end': 55181220,
                'start': 55181219,
                'type': 'SimpleInterval'
            },
            'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': '',
            'type': 'SequenceState'
        },
        'type': 'Allele'
    }),
    ("NC_000007.14:g.55181230_55181231insGGCT", {
        'location': {
            'interval': {
                'end': 55181230,
                'start': 55181230,
                'type': 'SimpleInterval'
            },
            'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'GGCT',
            'type': 'SequenceState'
        },
        'type': 'Allele'
    }),
    ("NC_000013.11:g.32331093_32331094dup", {
        'location': {
            'interval': {
                'end': 32331094,
                'start': 32331082,
                'type': 'SimpleInterval'
            },
            'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'TTTTTTTTTTTTTT',
            'type': 'SequenceState'
        },
        'type': 'Allele'
    }),
    ("NC_000013.11:g.32316467dup", {
        'location': {
            'interval': {
                'end': 32316467,
                'start': 32316466,
                'type': 'SimpleInterval'
            },
            'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'AA',
            'type': 'SequenceState'
        },
        'type': 'Allele'
    }),
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

# cases for testing to/from vcf records
vcf_variation_cases = [
    # multiple alleles
    (('1', '29881425', 'C', ['A', 'T']), {
        "type":
        "VariationSet",
        "members": [{
            '_id': 'ga4gh:VA.32Q-JABVQ7pecrW2RPlYYLxBeP6j6ZN_',
            'location': {
                'interval': {
                    'end': 29881425,
                    'start': 29881424,
                    'type': 'SimpleInterval'
                },
                'sequence_id': 'ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO',
                'type': 'SequenceLocation'
            },
            'state': {
                'sequence': 'A',
                'type': 'SequenceState'
            },
            'type': 'Allele'
        }, {
            '_id': 'ga4gh:VA.wUquYEB9ztjtsc9Xc4BnlpdQHOP4SvwN',
            'location': {
                'interval': {
                    'end': 29881425,
                    'start': 29881424,
                    'type': 'SimpleInterval'
                },
                'sequence_id': 'ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO',
                'type': 'SequenceLocation'
            },
            'state': {
                'sequence': 'T',
                'type': 'SequenceState'
            },
            'type': 'Allele'
        }]
    }),
    (('4', '147973044', 'C', ['A', 'T']), {
        "type":
        "VariationSet",
        "members": [{
            '_id': 'ga4gh:VA.I-Ol7SHx4gVYmvSDZ4tSrbKjQFQFnmAF',
            'location': {
                'interval': {
                    'end': 147973044,
                    'start': 147973043,
                    'type': 'SimpleInterval'
                },
                'sequence_id': 'ga4gh:SQ.HxuclGHh0XCDuF8x6yQrpHUBL7ZntAHc',
                'type': 'SequenceLocation'
            },
            'state': {
                'sequence': 'A',
                'type': 'SequenceState'
            },
            'type': 'Allele'
        }, {
            '_id': 'ga4gh:VA.SPrLYmcDQqlPOGYZWvE3AZ5skV9_EnXe',
            'location': {
                'interval': {
                    'end': 147973044,
                    'start': 147973043,
                    'type': 'SimpleInterval'
                },
                'sequence_id': 'ga4gh:SQ.HxuclGHh0XCDuF8x6yQrpHUBL7ZntAHc',
                'type': 'SequenceLocation'
            },
            'state': {
                'sequence': 'T',
                'type': 'SequenceState'
            },
            'type': 'Allele'
        }]
    }),
    (('7', '142149548', 'G', ['GT', 'GTT']), {
        "type":
        "VariationSet",
        "members": [{
            '_id': 'ga4gh:VA.y8-6HmQkyOLtgYCWozwfuSy8I8HbQocl',
            'location': {
                'interval': {
                    'end': 142149558,
                    'start': 142149548,
                    'type': 'SimpleInterval'
                },
                'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                'type': 'SequenceLocation'
            },
            'state': {
                'sequence': 'TTTTTTTTTTT',
                'type': 'SequenceState'
            },
            'type': 'Allele'
        }, {
            '_id': 'ga4gh:VA.z8z9fvWr1RHnOa_LjIr22rDnmXU7WkjZ',
            'location': {
                'interval': {
                    'end': 142149558,
                    'start': 142149548,
                    'type': 'SimpleInterval'
                },
                'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                'type': 'SequenceLocation'
            },
            'state': {
                'sequence': 'TTTTTTTTTTTT',
                'type': 'SequenceState'
            },
            'type': 'Allele'
        }]
    }),
    (('12', '49520727', 'A', ['AT', 'ATT', 'ATTT']), {
        "type":
        "VariationSet",
        "members": [{
            '_id': 'ga4gh:VA.y7R6DdlWtMVpe9H_1cKtuMS2egyf4inX',
            'location': {
                'interval': {
                    'end': 49520742,
                    'start': 49520727,
                    'type': 'SimpleInterval'
                },
                'sequence_id': 'ga4gh:SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl',
                'type': 'SequenceLocation'
            },
            'state': {
                'sequence': 'TTTTTTTTTTTTTTTT',
                'type': 'SequenceState'
            },
            'type': 'Allele'
        }, {
            '_id': 'ga4gh:VA.uZ_O9g5pteehH2E_gmjv4WkKK3lk4QOd',
            'location': {
                'interval': {
                    'end': 49520742,
                    'start': 49520727,
                    'type': 'SimpleInterval'
                },
                'sequence_id': 'ga4gh:SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl',
                'type': 'SequenceLocation'
            },
            'state': {
                'sequence': 'TTTTTTTTTTTTTTTTT',
                'type': 'SequenceState'
            },
            'type': 'Allele'
        }, {
            '_id': 'ga4gh:VA.kNYqY6TpWQ2EM400mW6fZrwRL1HZmz9S',
            'location': {
                'interval': {
                    'end': 49520742,
                    'start': 49520727,
                    'type': 'SimpleInterval'
                },
                'sequence_id': 'ga4gh:SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl',
                'type': 'SequenceLocation'
            },
            'state': {
                'sequence': 'TTTTTTTTTTTTTTTTTT',
                'type': 'SequenceState'
            },
            'type': 'Allele'
        }]
    }),
    # SNPs
    (('1', '92633', 'C', ['T']), {
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
        'state': {
            'type': 'SequenceState',
            'sequence': 'T'
        }
    }),
    (('1', '63002', 'A', ['G']), {
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
        'state': {
            'type': 'SequenceState',
            'sequence': 'G'
        }
    }),
    (('3', '166125806', 'G', ['A']), {
        '_id': 'ga4gh:VA.UAcSqA7iJtI36BysCOYXR8wPGVwUm4En',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX',
            'interval': {
                'type': 'SimpleInterval',
                'start': 166125805,
                'end': 166125806
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'A'
        }
    }),
    (('5', '1100700', 'G', ['A']), {
        '_id': 'ga4gh:VA.UrXpOOJGuCHAQgQ2sY85FsZ3CztIwI6C',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI',
            'interval': {
                'type': 'SimpleInterval',
                'start': 1100699,
                'end': 1100700
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'A'
        }
    }),
    (('7', '108934752', 'T', ['C']), {
        '_id': 'ga4gh:VA.pwoWgFDuLdigDfQRTEU7gCdSB4DoHED1',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
            'interval': {
                'type': 'SimpleInterval',
                'start': 108934751,
                'end': 108934752
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'C'
        }
    }),
    (('13', '59140800', 'G', ['A']), {
        '_id': 'ga4gh:VA.Ii5EzfJNVealFy9Wc7UMjd4WrzqSpbmg',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
            'interval': {
                'type': 'SimpleInterval',
                'start': 59140799,
                'end': 59140800
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'A'
        }
    }),
    (('X', '127226377', 'A', ['G']), {
        '_id': 'ga4gh:VA.nyaA37BGIVLcP6LPfpT1GKg-nfLhxPOD',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP',
            'interval': {
                'type': 'SimpleInterval',
                'start': 127226376,
                'end': 127226377
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'G'
        }
    }),
    # ins
    (('1', '72297', 'G', ['GTAT']), {
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
        'state': {
            'type': 'SequenceState',
            'sequence': 'TATTATT'
        }
    }),
    (('3', '166114496', 'T', ['TA']), {
        '_id': 'ga4gh:VA.LVqeAmYsZxO3kYb-lvmYOrqnoY1NbKWd',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX',
            'interval': {
                'type': 'SimpleInterval',
                'start': 166114496,
                'end': 166114496
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'A'
        }
    }),
    (('5', '127031196', 'C', ['CG']), {
        '_id': 'ga4gh:VA.4lIsTJ1zxjQ7UjXTSnCLb-HF2jw733NK',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI',
            'interval': {
                'type': 'SimpleInterval',
                'start': 127031196,
                'end': 127031201
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'GGGGGG'
        }
    }),
    (('7', '156408692', 'C', ['CAT']), {
        '_id': 'ga4gh:VA.uvVtM0pr6hFMLfKTurDcdFj_3rRHbBqp',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
            'interval': {
                'type': 'SimpleInterval',
                'start': 156408692,
                'end': 156408703
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'ATATATATATATA'
        }
    }),
    (('10', '106329140', 'C', ['CATTT']), {
        '_id': 'ga4gh:VA.7tZ-MNUx-ZXdOElypDtzWJ1jL5JbMYgE',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.ss8r_wB0-b9r44TQTMmVTI92884QvBiB',
            'interval': {
                'type': 'SimpleInterval',
                'start': 106329140,
                'end': 106329143
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'ATTTATT'
        }
    }),
    (('Y', '22304601', 'G', ['GA']), {
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
        'state': {
            'type': 'SequenceState',
            'sequence': 'AAAAAAA'
        }
    }),
    # del
    (('4', '116619313', 'GT', ['G']), {
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
        'state': {
            'type': 'SequenceState',
            'sequence': ''
        }
    }),
    (('5', '81229982', 'TG', ['T']), {
        '_id': 'ga4gh:VA.UU86eospaxRFVjtkX0VHKBbcNenwDfU3',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI',
            'interval': {
                'type': 'SimpleInterval',
                'start': 81229982,
                'end': 81229983
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': ''
        }
    }),
    (('5', '81230677', 'AG', ['A']), {
        '_id': 'ga4gh:VA.xJctu09E74WfvRu-TvtjesE8S_spoKoB',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI',
            'interval': {
                'type': 'SimpleInterval',
                'start': 81230677,
                'end': 81230684
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'GGGGGG'
        }
    }),
    (('5', '81233450', 'TGGA', ['T']), {
        '_id': 'ga4gh:VA.JxnxekfJP329pXfBum6U8a5u0NSNFmnp',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI',
            'interval': {
                'type': 'SimpleInterval',
                'start': 81233450,
                'end': 81233457
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'GGAG'
        }
    }),
    (('5', '81233494', 'GCTGA', ['G']), {
        '_id': 'ga4gh:VA.-lNWDXnYR4hs2u6fFncH1ZOqILu4xZT8',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI',
            'interval': {
                'type': 'SimpleInterval',
                'start': 81233494,
                'end': 81233501
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'CTG'
        }
    }),
    (('8', '62540756', 'TCTC', ['T']), {
        '_id': 'ga4gh:VA.oq5H3D7ea1BLpVONa-UmLhDiKug7kOaU',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.209Z7zJ-mFypBEWLk4rNC6S_OxY5p7bs',
            'interval': {
                'type': 'SimpleInterval',
                'start': 62540756,
                'end': 62540761
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'CT'
        }
    }),
    (('10', '106332351', 'CAT', ['C']), {
        '_id': 'ga4gh:VA.xti_mlNXAekytcegCcaqk9csAQ5f5pYw',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.ss8r_wB0-b9r44TQTMmVTI92884QvBiB',
            'interval': {
                'type': 'SimpleInterval',
                'start': 106332351,
                'end': 106332355
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'AT'
        }
    }),
    (('17', '29204173', 'TACA', ['T']), {
        '_id': 'ga4gh:VA.mId9FgDqwHkP3mBt1lgfSr2-bhCJaxGg',
        'type': 'Allele',
        'location': {
            'type': 'SequenceLocation',
            'sequence_id': 'ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7',
            'interval': {
                'type': 'SimpleInterval',
                'start': 29204173,
                'end': 29204179
            }
        },
        'state': {
            'type': 'SequenceState',
            'sequence': 'ACA'
        }
    })
]


@pytest.mark.parametrize("record,expected", vcf_variation_cases)
def test_from_vcf_record(tlr_norm, record, expected):
    """Test Translator._from_vcf_record"""
    variation_from_vcf = tlr_norm._from_vcf_record(record[0], record[1], record[2], record[3], assembly_name='GRCh38')
    assert variation_from_vcf.as_dict() == expected


@pytest.fixture(scope="session")
def vcf_file_in():
    """Provide input file fixture to test Translator._from_vcf"""
    file_string = '##fileformat=VCFv4.3\n' + \
                  '##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta\n' + \
                  '##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,' + \
                  'species="Homo sapiens",taxonomy=x>\n' + \
                  '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n' + \
                  '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n' + \
                  '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">\n' + \
                  '##FILTER=<ID=q10,Description="Quality below 10">\n' + \
                  '##FILTER=<ID=s50,Description="Less than 50% of samples have data">\n' + \
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' + \
                  '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n' + \
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003\n' + \
                  '20\t14370\trs6054257\tG\tA\t29\tPASS\tDP=14;AF=0.5;DB\tGT:DP\t0/0:1\t0/1:8\t1/1:5\n'
    file_bytes = file_string.encode('utf-8')
    fp = TemporaryFile()
    fp.write(file_bytes)
    fp.seek(0)
    return fp


@pytest.fixture(scope="session")
def vcf_from_file_expected():
    """Provide expected output to test Translator._from_vcf"""
    return {
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
        'state': {
            'type': 'SequenceState',
            'sequence': 'A'
        }
    },


def test_from_vcf_file(tlr_norm, vcf_file_in, vcf_from_file_expected):
    """Test Translator._from_vcf"""
    alleles = tlr_norm._from_vcf(vcf_file_in)
    for i in range(len(alleles)):
        assert alleles[i].as_dict() == vcf_from_file_expected[i]


@pytest.fixture(scope="session")
def vcf_to_file():
    return vcf_variation_cases


def test_to_vcf_file(tlr_norm, vcf_to_file):
    """Test Translator._to_vcf"""
    # alleles = [tlr_norm._from_vrs(i) for j in vcf_to_file for i in j[1]]
    alleles = [tlr_norm._from_vrs(i[1]) for i in vcf_to_file]
    outfile_path = "test_out.vcf"
    tlr_norm._to_vcf(alleles, outfile_path)

    # get actual rows
    with open(outfile_path) as f:
        outfile_lines = list(f.readlines())

    # build expected rows -- sort to match output rules
    def format_as_vcf_row(tup):
        return f'{tup[0]}\t{tup[1]}\t.\t{tup[2]}\t{",".join(tup[3])}\t.\t.\t.\n'
    expected = [i[0] for i in vcf_to_file]
    expected.sort(key=lambda r: (int(r[0]), int(r[1])) if r[0] not in ('X', 'Y') else (ord(r[0]), int(r[1])))

    for i in range(len(vcf_to_file)):
        assert outfile_lines[i + 15] == format_as_vcf_row(expected[i])

    remove(outfile_path)
