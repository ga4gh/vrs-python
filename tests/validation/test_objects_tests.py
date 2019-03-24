import os

import pytest
import yaml

import vmc


validation_fn = os.path.join(os.path.dirname(__file__), "..", "..", "vmc-spec", "tests", "objects.yaml")
validation_tests = yaml.load(open(validation_fn), Loader=yaml.SafeLoader)


@pytest.mark.parametrize("test", validation_tests["Interval"])
def test_Interval(test):
    o = vmc.models.Interval(**test["in"])
    assert test["out"]["serialize"] == vmc.serialize(o)
    assert test["out"]["computed_id"] == vmc.computed_id(o)


@pytest.mark.parametrize("test", validation_tests["Location"])
def test_Interval(test):
    o = vmc.models.Location(**test["in"])
    assert test["out"]["serialize"] == vmc.serialize(o)
    assert test["out"]["computed_id"] == vmc.computed_id(o)


@pytest.mark.parametrize("test", validation_tests["Allele"])
def test_Allele(test):
    o = vmc.models.Allele(**test["in"])
    assert test["out"]["serialize"] == vmc.serialize(o)
    assert test["out"]["computed_id"] == vmc.computed_id(o)


@pytest.mark.parametrize("test", validation_tests["Haplotype"])
def test_Haplotype(test):
    o = vmc.models.Haplotype(**test["in"])
    assert test["out"]["serialize"] == vmc.serialize(o)
    assert test["out"]["computed_id"] == vmc.computed_id(o)


@pytest.mark.parametrize("test", validation_tests["Genotype"])
def test_Genotype(test):
    o = vmc.models.Genotype(**test["in"])
    assert test["out"]["serialize"] == vmc.serialize(o)
    assert test["out"]["computed_id"] == vmc.computed_id(o)
