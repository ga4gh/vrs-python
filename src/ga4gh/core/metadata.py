"""Provide shared metadata types for GA4GH GKS models."""

from enum import Enum
from typing import Any, ClassVar

from pydantic.json_schema import GenerateJsonSchema, JsonSchemaMode


class Maturity(str, Enum):
    """Maturity levels for GA4GH product features."""

    DRAFT = "draft"
    TRIAL_USE = "trial use"
    NORMATIVE = "normative"
    DEPRECATED = "deprecated"


class GKSMaturityMixin:
    """Provide maturity metadata for a GA4GH GKS model."""

    _maturity: ClassVar[Maturity]

    @classmethod
    def maturity(cls) -> Maturity:
        """Return the GKS maturity level for the model."""
        return cls._maturity


class GKSSchemaMixin:
    """Provide a canonical JSON Schema identifier for a GA4GH GKS model."""

    _schema_base_uri: ClassVar[str] = "https://w3id.org/ga4gh/schema"
    _product_name: ClassVar[str]
    _product_version: ClassVar[str]

    @classmethod
    def schema_id(cls) -> str:
        """Return the canonical JSON Schema identifier for the model."""
        return f"{cls._schema_base_uri}/{cls._product_name}/{cls._product_version}/json/{cls.__name__}"


class GKSMetadataMixin(GKSMaturityMixin, GKSSchemaMixin):
    """Provide maturity and schema metadata for a concrete GKS model."""

    @classmethod
    def model_json_schema(
        cls,
        by_alias: bool = True,
        ref_template: str = "#/$defs/{model}",
        schema_generator: type[GenerateJsonSchema] = GenerateJsonSchema,
        mode: JsonSchemaMode = "validation",
    ) -> dict[str, Any]:
        """Generate JSON Schema with GKS metadata."""
        schema = super().model_json_schema(
            by_alias=by_alias,
            ref_template=ref_template,
            schema_generator=schema_generator,
            mode=mode,
        )

        schema["$id"] = cls.schema_id()
        schema["maturity"] = cls.maturity().value

        ga4gh_class = getattr(cls, "ga4gh", None)
        if not ga4gh_class:
            return schema

        ga4gh_metadata = {}

        if prefix := getattr(ga4gh_class, "prefix", None):
            ga4gh_metadata["prefix"] = prefix

        if inherent := getattr(ga4gh_class, "inherent", None):
            ga4gh_metadata["inherent"] = list(inherent)

        if ga4gh_metadata:
            schema["ga4gh"] = ga4gh_metadata

        return schema
