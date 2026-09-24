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
        / "gkm-core"
        / "schema"
        / "gkm-core"
        / "gkm-core-source.yaml",
        SUBMODULES_DIR / "submodules" / "gkm-core" / "schema" / "gkm-core" / "json",
    ),
    (
        vrs_models,
        SUBMODULES_DIR / "schema" / "vrs" / "vrs-source.yaml",
        SUBMODULES_DIR / "schema" / "vrs" / "json",
    ),
)


def _concrete_model_params():
    """Return concrete model metadata discovered from JSON Schema files.

    :returns: Pytest parameters for concrete GKS models.
    """
    params = []
    for model_module, _, json_dir in SCHEMAS:
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
    """Return abstract model metadata found only in source schemas.

    :returns: Pytest parameters for abstract GKS models and source definitions.
    """
    params = []
    for model_module, source_path, _ in SCHEMAS:
        schema_params = []
        with source_path.open() as source_file:
            definitions = yaml.safe_load(source_file)["$defs"]

        for name, definition in definitions.items():
            if definition.get("abstract") is True:
                schema_params.append(
                    pytest.param(getattr(model_module, name), definition, id=name)
                )

        assert schema_params, f"No abstract models discovered in {source_path}"

        params.extend(schema_params)
    return params


def _abstract_schema_model_params():
    """Return abstract model metadata discovered from source schemas.

    :returns: Pytest parameters for abstract GKS models.
    """
    params = []
    for model_module, source_path, _ in SCHEMAS:
        with source_path.open() as source_file:
            definitions = yaml.safe_load(source_file)["$defs"]
        for name, definition in definitions.items():
            if definition.get("abstract") is True:
                params.append(pytest.param(getattr(model_module, name), id=name))
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
    """Verify abstract models expose their source-defined maturity.

    :param model: Abstract Pydantic model.
    :param definition: Corresponding source schema definition.
    """
    assert "_maturity" in model.__dict__
    assert model.maturity() == Maturity(definition["maturity"])


@pytest.mark.parametrize("model", _abstract_schema_model_params())
def test_abstract_model_schema_metadata(model):
    """Verify abstract models emit the abstract schema keyword.

    :param model: Abstract Pydantic model.
    """
    assert model.model_json_schema()["abstract"] is True


@pytest.mark.parametrize("model", _abstract_schema_model_params())
def test_abstract_models_cannot_be_instantiated(model):
    """Verify abstract models reject direct construction.

    :param model: Abstract Pydantic model.
    """
    kwargs = {} if model is core_models.Element else {"type": "test"}
    with pytest.raises(ValueError, match="abstract and cannot be instantiated"):
        model(**kwargs)
