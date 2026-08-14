"""Test model metadata against the GKS source and JSON schemas."""

import json
from pathlib import Path

import pytest
import yaml

from ga4gh.core import core_models
from ga4gh.core.metadata import Maturity
from ga4gh.vrs import models as vrs_models

SUBMODULES_DIR = Path(__file__).parents[2] / "submodules" / "vrs"
SCHEMAS = (
    (
        core_models,
        SUBMODULES_DIR
        / "submodules"
        / "gks-core"
        / "schema"
        / "gks-core"
        / "gks-core-source.yaml",
        SUBMODULES_DIR / "submodules" / "gks-core" / "schema" / "gks-core" / "json",
    ),
    (
        vrs_models,
        SUBMODULES_DIR / "schema" / "vrs" / "vrs-source.yaml",
        SUBMODULES_DIR / "schema" / "vrs" / "json",
    ),
)


def _concrete_model_params():
    """Return concrete model metadata discovered from JSON Schema files."""
    params = []
    for model_module, _, json_dir in SCHEMAS:
        schema_params = []
        for schema_path in sorted(json_dir.iterdir()):
            model = getattr(model_module, schema_path.name, None)
            if model is None:
                continue  # date and datetime use standard-library classes
            with schema_path.open() as schema_file:
                schema = json.load(schema_file)
            schema_params.append(pytest.param(model, schema, id=schema["title"]))
        assert schema_params, f"No concrete models discovered in {json_dir}"
        params.extend(schema_params)
    return params


def _abstract_model_params():
    """Return abstract model metadata found only in source schemas."""
    params = []
    for model_module, source_path, json_dir in SCHEMAS:
        schema_params = []
        with source_path.open() as source_file:
            definitions = yaml.safe_load(source_file)["$defs"]
        concrete_names = {path.name for path in json_dir.iterdir()}
        for name, definition in definitions.items():
            if name not in concrete_names and "heritableProperties" in definition:
                schema_params.append(
                    pytest.param(getattr(model_module, name), definition, id=name)
                )
        assert schema_params, f"No abstract models discovered in {source_path}"
        params.extend(schema_params)
    return params


@pytest.mark.parametrize(("model", "schema"), _concrete_model_params())
def test_concrete_model_metadata(model, schema):
    """Concrete model metadata matches its generated JSON Schema."""
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
    """Abstract models expose source-defined maturity but no schema identifier."""
    assert model.maturity() == Maturity(definition["maturity"])
    assert not hasattr(model, "schema_id")
