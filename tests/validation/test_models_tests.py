import os

import pytest
import yaml

import vmc
from vmc.serialize import serialize_cj, serialize_vmc


validation_fn = os.path.join(os.path.dirname(__file__), "..", "..", "vmc-spec", "tests", "models.yaml")
validation_tests = yaml.load(open(validation_fn), Loader=yaml.SafeLoader)


@pytest.mark.parametrize("test", validation_tests["SimpleRegion"])
def test_SimpleRegion(test):
    o = vmc.models.SimpleRegion(**test["in"])
    assert test["out"]["serialize_cj"] == serialize_cj(o).decode("utf-8")
    assert test["out"]["serialize_vmc"] == serialize_vmc(o)


# @pytest.mark.parametrize("test", validation_tests["Location"])
# def test_Interval(test):
#     o = vmc.models.Location(**test["in"])
#     assert test["out"]["serialize"] == vmc.serialize(o)
#     assert test["out"]["computed_id"] == vmc.computed_id(o)
# 
# 
# @pytest.mark.parametrize("test", validation_tests["Allele"])
# def test_Allele(test):
#     o = vmc.models.Allele(**test["in"])
#     assert test["out"]["serialize"] == vmc.serialize(o)
#     assert test["out"]["computed_id"] == vmc.computed_id(o)
# 
# 
# @pytest.mark.parametrize("test", validation_tests["Haplotype"])
# def test_Haplotype(test):
#     o = vmc.models.Haplotype(**test["in"])
#     assert test["out"]["serialize"] == vmc.serialize(o)
#     assert test["out"]["computed_id"] == vmc.computed_id(o)
# 
# 
# @pytest.mark.parametrize("test", validation_tests["Genotype"])
# def test_Genotype(test):
#     o = vmc.models.Genotype(**test["in"])
#     assert test["out"]["serialize"] == vmc.serialize(o)
#     assert test["out"]["computed_id"] == vmc.computed_id(o)
