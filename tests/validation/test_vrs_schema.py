"""test that VRS Python model structures match VRS Schema
"""
import yaml
from pathlib import Path
from ga4gh.vrs import models

ROOT_DIR = Path(__file__).parents[2]
VRS_SCHEMA_DIR = ROOT_DIR / 'submodules' / 'vrs' / 'schema'

with open(VRS_SCHEMA_DIR / 'vrs.yaml') as vrs_yaml:
    VRS_SCHEMA = yaml.safe_load(vrs_yaml)

VRS_CONCRETE_CLASSES = set()
VRS_PRIMITIVES = set()

for vrs_class in VRS_SCHEMA['$defs']:
    cls_def = VRS_SCHEMA['$defs'][vrs_class]
    if 'properties' in cls_def:
        VRS_CONCRETE_CLASSES.add(vrs_class)
    elif cls_def.get('type') in ['array', 'int', 'str']:
        VRS_PRIMITIVES.add(vrs_class)

NOT_IMPLEMENTED = ['Adjacency', 'Haplotype']  # Use this to skip testing of not-implemented classes
                                              # TODO: Remove this once 2.0 models at beta


def test_schema_models_exist():
    """test that VRS Python covers the models defined by VRS
    """
    for vrs_class in VRS_CONCRETE_CLASSES | VRS_PRIMITIVES:
        if vrs_class in NOT_IMPLEMENTED:
            continue
        assert getattr(models, vrs_class, False)


def test_schema_class_fields_are_valid():
    """test that VRS Python model fields match the VRS specification
    """
    for vrs_class in VRS_CONCRETE_CLASSES:
        if vrs_class in NOT_IMPLEMENTED:
            continue
        schema_fields = set(VRS_SCHEMA['$defs'][vrs_class]['properties'])
        pydantic_model = getattr(models, vrs_class)
        assert set(pydantic_model.__fields__) == schema_fields, vrs_class


def test_model_keys_are_valid():
    """test that digest keys on Value Objects are valid and sorted
    """
    for vrs_class in VRS_CONCRETE_CLASSES:
        if vrs_class in NOT_IMPLEMENTED:
            continue
        if VRS_SCHEMA['$defs'][vrs_class].get('ga4ghDigest', {}).get('keys', None) is None:
            continue
        pydantic_model = getattr(models, vrs_class)
        try:
            pydantic_model_digest_keys = pydantic_model.ga4gh.keys
        except AttributeError:
            raise AttributeError(vrs_class)
        assert set(pydantic_model_digest_keys) == set(VRS_SCHEMA['$defs'][vrs_class]['ga4ghDigest']['keys']), vrs_class
        assert pydantic_model_digest_keys == sorted(pydantic_model.ga4gh.keys), vrs_class
