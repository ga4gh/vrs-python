import pytest

from ga4gh.vrs import models
from ga4gh.vrs.extras.variation_normalizer_rest_dp import VariationNormalizerRESTDataProxy
from tests.extras.test_translator import hgvs_tests


@pytest.fixture(scope="module")
def variation_norm_rest_dp():
    return VariationNormalizerRESTDataProxy()


@pytest.mark.parametrize("expected,vo_as_dict", hgvs_tests)
@pytest.mark.vcr
def test_rest_dp_to_hgvs(variation_norm_rest_dp, expected, vo_as_dict):
    vo = models.Allele(**vo_as_dict)
    resp = variation_norm_rest_dp.to_hgvs(vo)
    assert resp == [expected]
