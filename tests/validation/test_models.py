import os

import pytest
import yaml

from ga4gh.vr import models, ga4gh_serialize


validation_fn = os.path.join(os.path.dirname(__file__), "data", "models.yaml")
validation_tests = yaml.load(open(validation_fn), Loader=yaml.SafeLoader)


@pytest.mark.parametrize("test", validation_tests["SimpleInterval"])
def test_SimpleInterval(test):
    o = models.SimpleInterval(**test["in"])
    assert test["out"]["serialize"] == ga4gh_serialize(o).decode()

@pytest.mark.parametrize("test", validation_tests["NestedInterval"])
def test_NestedInterval(test):
    o = models.NestedInterval(**test["in"])
    assert test["out"]["serialize"] == ga4gh_serialize(o).decode()

@pytest.mark.parametrize("test", validation_tests["SequenceLocation"])
def test_SequenceLocation(test):
    o = models.SequenceLocation(**test["in"])
    assert test["out"]["serialize"] == ga4gh_serialize(o).decode()

@pytest.mark.parametrize("test", validation_tests["Allele"])
def test_Allele(test):
    o = models.Allele(**test["in"])
    assert test["out"]["serialize"] == ga4gh_serialize(o).decode()
