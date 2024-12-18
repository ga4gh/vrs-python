"""Test that VRS-Python Pydantic models match VRS and GKS-Common schemas"""

from enum import Enum
import json
from pathlib import Path

import pytest
from pydantic import BaseModel

from ga4gh.core import core_models
from ga4gh.vrs import models as vrs_models


class GKSSchema(str, Enum):
    """Enum for GKS schema"""

    VRS = "vrs"
    CORE = "core"


class GKSSchemaMapping(BaseModel):
    """Model for representing GKS Schema concrete classes, primitives, and schema"""

    base_classes: set = set()
    concrete_classes: set = set()
    primitives: set = set()
    schema: dict = {}


def _update_gks_schema_mapping(
    f_path: Path, gks_schema_mapping: GKSSchemaMapping
) -> None:
    """Update ``gks_schema_mapping`` properties

    :param f_path: Path to JSON Schema file
    :param gks_schema_mapping: GKS schema mapping to update
    """
    with f_path.open() as rf:
        cls_def = json.load(rf)

    spec_class = cls_def["title"]
    gks_schema_mapping.schema[spec_class] = cls_def

    if "properties" in cls_def:
        gks_schema_mapping.concrete_classes.add(spec_class)
    elif cls_def.get("type") in {"array", "integer", "string"}:
        gks_schema_mapping.primitives.add(spec_class)
    else:
        gks_schema_mapping.base_classes.add(spec_class)


GKS_SCHEMA_MAPPING = {gks: GKSSchemaMapping() for gks in GKSSchema}
SUBMODULES_DIR = Path(__file__).parents[2] / "submodules" / "vrs"


# Get vrs classes
vrs_mapping = GKS_SCHEMA_MAPPING[GKSSchema.VRS]
for f in (SUBMODULES_DIR / "schema" / "vrs" / "json").glob("*"):
    _update_gks_schema_mapping(f, vrs_mapping)

# Get core classes
core_mapping = GKS_SCHEMA_MAPPING[GKSSchema.CORE]
for f in (SUBMODULES_DIR / "submodules" / "gks-core" / "schema" / "gks-core" / "json").glob("*"):
    _update_gks_schema_mapping(f, core_mapping)


@pytest.mark.parametrize(
    ("gks_schema", "pydantic_models"),
    [
        (GKSSchema.VRS, vrs_models),
        (GKSSchema.CORE, core_models),
    ],
)
def test_schema_models_in_pydantic(gks_schema, pydantic_models):
    """Ensure that each schema model has corresponding Pydantic model"""
    mapping = GKS_SCHEMA_MAPPING[gks_schema]
    for schema_model in (
        mapping.base_classes | mapping.concrete_classes | mapping.primitives
    ):
        assert getattr(pydantic_models, schema_model, False), schema_model


@pytest.mark.parametrize(
    ("gks_schema", "pydantic_models"),
    [
        (GKSSchema.VRS, vrs_models),
        (GKSSchema.CORE, core_models),
    ],
)
def test_schema_class_fields(gks_schema, pydantic_models):
    """Check that each schema model properties exist and are required in corresponding
    Pydantic model, and validate required properties
    """
    mapping = GKS_SCHEMA_MAPPING[gks_schema]
    for schema_model in mapping.concrete_classes:
        schema_properties = mapping.schema[schema_model]["properties"]
        pydantic_model = getattr(pydantic_models, schema_model)
        assert set(pydantic_model.model_fields) == set(schema_properties), schema_model

        required_schema_fields = set(mapping.schema[schema_model]["required"])

        for prop, property_def in schema_properties.items():
            pydantic_model_field_info = pydantic_model.model_fields[prop]
            pydantic_field_required = pydantic_model_field_info.is_required()

            if prop in required_schema_fields:
                if prop != "type":
                    assert pydantic_field_required, f"{pydantic_model}.{prop}"
            else:
                assert not pydantic_field_required, f"{pydantic_model}.{prop}"

            if "description" in property_def:
                if prop != "code":  # special exception
                    assert property_def["description"].replace("'", "\"") == pydantic_model_field_info.description.replace("'", "\""), f"{pydantic_model}.{prop}"
            else:
                assert pydantic_model_field_info.description is None, f"{pydantic_model}.{prop}"


def test_ga4gh_keys():
    """Ensure ga4gh inherit defined in schema model exist in corresponding Pydantic model"""
    vrs_mapping = GKS_SCHEMA_MAPPING[GKSSchema.VRS]
    for vrs_class in vrs_mapping.concrete_classes:
        if (
            vrs_mapping.schema[vrs_class].get("ga4gh", {}).get("inherit", None)
            is None
        ):
            continue

        pydantic_model = getattr(vrs_models, vrs_class)

        try:
            pydantic_model_digest_keys = pydantic_model.ga4gh.keys
        except AttributeError as e:
            raise AttributeError(vrs_class) from e

        assert set(pydantic_model_digest_keys) == set(
            vrs_mapping.schema[vrs_class]["ga4gh"]["inherit"]
        ), vrs_class
        assert pydantic_model_digest_keys == sorted(
            pydantic_model.ga4gh.keys
        ), vrs_class
