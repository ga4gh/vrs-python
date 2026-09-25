"""Test model metadata against the GKM JSON schemas."""

import json
from pathlib import Path

import pytest
from pydantic import RootModel

from ga4gh.core import core_models
from ga4gh.core.metadata import (
    GKMMaturityMixin,
    GKMMetadataMixin,
    GKMSchemaMixin,
    GKSMaturityMixin,
    GKSMetadataMixin,
    GKSSchemaMixin,
    Maturity,
)
from ga4gh.vrs import models as vrs_models

SUBMODULES_DIR = Path(__file__).parents[2] / "submodules" / "vrs"
SCHEMAS = (
    (
        core_models,
        SUBMODULES_DIR / "submodules" / "gkm-core" / "schema" / "gkm-core" / "json",
    ),
    (
        vrs_models,
        SUBMODULES_DIR / "schema" / "vrs" / "json",
    ),
)


@pytest.mark.parametrize(
    ("deprecated_model", "canonical_model"),
    [
        (GKSMaturityMixin, GKMMaturityMixin),
        (GKSSchemaMixin, GKMSchemaMixin),
        (GKSMetadataMixin, GKMMetadataMixin),
        (core_models.GKSCoreMetadataMixin, core_models.GKMCoreMetadataMixin),
    ],
)
def test_gks_models_are_deprecated(deprecated_model, canonical_model):
    """GKS model names remain available as deprecated aliases."""
    with pytest.deprecated_call():
        deprecated_model()
    assert issubclass(deprecated_model, canonical_model)


def _concrete_model_params():
    """Return concrete model metadata discovered from JSON Schema files.

    :returns: Pytest parameters for concrete GKM models.
    """
    params = []
    for model_module, json_dir in SCHEMAS:
        schema_params = []
        for schema_path in sorted(json_dir.iterdir()):
            model = getattr(model_module, schema_path.name, None)
            if model is None:
                continue  # date and datetime use standard-library classes

            with schema_path.open() as schema_file:
                schema = json.load(schema_file)

            if schema.get("abstract") is True:
                continue

            schema_params.append(pytest.param(model, schema, id=schema["title"]))
        assert schema_params, f"No concrete models discovered in {json_dir}"
        params.extend(schema_params)
    return params


def _abstract_model_params():
    """Return abstract model metadata from JSON Schema files.

    :returns: Pytest parameters for abstract GKM models and JSON definitions.
    """
    params = []
    for model_module, json_dir in SCHEMAS:
        schema_params = []
        for schema_path in sorted(json_dir.iterdir()):
            with schema_path.open() as schema_file:
                definition = json.load(schema_file)
            if definition.get("abstract") is True:
                schema_params.append(
                    pytest.param(
                        getattr(model_module, schema_path.name),
                        definition,
                        id=schema_path.name,
                    )
                )

        assert schema_params, f"No abstract models discovered in {json_dir}"

        params.extend(schema_params)
    return params


@pytest.mark.parametrize(("model", "schema"), _concrete_model_params())
def test_concrete_model_metadata(model, schema):
    """Verify concrete model metadata matches generated JSON Schema.

    :param model: Concrete Pydantic model.
    :param schema: Corresponding generated JSON Schema.
    """
    assert model.schema_id() == schema["$id"]
    assert model.maturity() == Maturity(schema["maturity"])
    generated_schema = model.model_json_schema()
    assert generated_schema["$id"] == schema["$id"]
    assert generated_schema["maturity"] == schema["maturity"]
    if ga4gh_metadata := schema.get("ga4gh"):
        assert generated_schema["ga4gh"].get("prefix") == ga4gh_metadata.get("prefix")
        assert set(generated_schema["ga4gh"]["inherent"]) == set(
            ga4gh_metadata["inherent"]
        )
    else:
        assert "ga4gh" not in generated_schema


@pytest.mark.parametrize(("model", "definition"), _abstract_model_params())
def test_abstract_model_metadata(model, definition):
    """Verify abstract models expose JSON Schema metadata.

    :param model: Abstract Pydantic model.
    :param definition: Corresponding JSON Schema definition.
    """
    assert "_maturity" in model.__dict__
    assert model.maturity() == Maturity(definition["maturity"])
    generated_schema = model.model_json_schema()
    assert generated_schema["$id"] == definition["$id"]
    assert generated_schema["maturity"] == definition["maturity"]
    assert generated_schema["abstract"] is True
    if issubclass(model, RootModel):
        # These are public compatibility adapters for the former sealed unions.
        # Pydantic adds a discriminator mapping and local $defs references, whereas
        # the published abstract schemas use portable references.
        assert generated_schema["discriminator"]["propertyName"] == "type"
        assert len(generated_schema["oneOf"]) == len(definition["oneOf"])
    else:
        assert generated_schema.get("discriminator") == definition.get("discriminator")
        assert generated_schema.get("oneOf") == definition.get("oneOf")


@pytest.mark.parametrize(
    ("model", "member", "payload"),
    [
        (
            vrs_models.Variation,
            vrs_models.CopyNumberChange,
            {
                "type": "CopyNumberChange",
                "location": "ga4gh:VSL.test",
                "copyChange": "loss",
            },
        ),
        (
            vrs_models.MolecularVariation,
            vrs_models.Allele,
            {
                "type": "Allele",
                "location": "ga4gh:VSL.test",
                "state": {"type": "LiteralSequenceExpression", "sequence": "A"},
            },
        ),
        (
            vrs_models.SystemicVariation,
            vrs_models.CopyNumberCount,
            {"type": "CopyNumberCount", "location": "ga4gh:VSL.test", "copies": 2},
        ),
        (
            vrs_models.SequenceExpression,
            vrs_models.LiteralSequenceExpression,
            {"type": "LiteralSequenceExpression", "sequence": "A"},
        ),
        (
            vrs_models.Location,
            vrs_models.SequenceLocation,
            {
                "type": "SequenceLocation",
                "sequenceReference": "SQ.test",
                "start": 1,
                "end": 2,
            },
        ),
    ],
)
def test_abstract_vrs_models_dispatch_typed_payloads(model, member, payload):
    """Abstract VRS models dispatch typed payloads to their concrete members."""
    result = model.model_validate(payload)
    assert isinstance(result.root, member)
    assert isinstance(model(root=payload).root, member)
