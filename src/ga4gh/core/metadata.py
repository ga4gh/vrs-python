"""Provide shared metadata types for GA4GH GKM models."""

from enum import Enum
from typing import Any, ClassVar

from pydantic.json_schema import GenerateJsonSchema, JsonSchemaMode
from typing_extensions import deprecated


class Maturity(str, Enum):
    """Maturity levels for GA4GH product features."""

    DRAFT = "draft"
    TRIAL_USE = "trial use"
    NORMATIVE = "normative"
    DEPRECATED = "deprecated"


class GKMMaturityMixin:
    """Provide maturity metadata for a GA4GH GKM model."""

    _maturity: ClassVar[Maturity]

    @classmethod
    def maturity(cls) -> Maturity:
        """Return the GKM maturity level for the model."""
        return cls._maturity


class GKMSchemaMixin:
    """Provide a canonical JSON Schema identifier for a GA4GH GKM model."""

    _schema_base_uri: ClassVar[str] = "https://w3id.org/ga4gh/schema"
    _product_name: ClassVar[str]
    _product_version: ClassVar[str]

    @classmethod
    def schema_id(cls) -> str:
        """Return the canonical JSON Schema identifier for the model."""
        return f"{cls._schema_base_uri}/{cls._product_name}/{cls._product_version}/json/{cls.__name__}"


class GKMMetadataMixin(GKMMaturityMixin, GKMSchemaMixin):
    """Provide maturity and schema metadata for a GKM model."""

    _abstract: ClassVar[bool] = False

    @staticmethod
    def apply_schema_metadata(
        model_class: type, schema: dict[str, Any]
    ) -> dict[str, Any]:
        """Add GKM metadata to a generated JSON Schema.

        :param model_class: Pydantic model class that produced the schema.
        :param schema: Generated JSON Schema to annotate.
        :returns: The annotated JSON Schema.
        """
        schema["$id"] = model_class.schema_id()
        schema["maturity"] = model_class.maturity().value

        if model_class.__dict__.get("_abstract", False):
            schema["abstract"] = True

        # GA4GH identifier metadata is optional and applies only when declared.
        ga4gh_class = getattr(model_class, "ga4gh", None)
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

    @classmethod
    def model_json_schema(
        cls,
        by_alias: bool = True,
        ref_template: str = "#/$defs/{model}",
        schema_generator: type[GenerateJsonSchema] = GenerateJsonSchema,
        mode: JsonSchemaMode = "validation",
    ) -> dict[str, Any]:
        """Generate JSON Schema with GKM metadata.

        :param by_alias: Whether to use field aliases.
        :param ref_template: Template for schema references.
        :param schema_generator: Pydantic schema generator class.
        :param mode: Pydantic schema generation mode.
        :returns: JSON Schema annotated with GKM metadata.
        """
        schema = super().model_json_schema(
            by_alias=by_alias,
            ref_template=ref_template,
            schema_generator=schema_generator,
            mode=mode,
        )

        return cls.apply_schema_metadata(cls, schema)


@deprecated("GKSMaturityMixin is deprecated; use GKMMaturityMixin instead.")
class GKSMaturityMixin(GKMMaturityMixin):
    """Deprecated alias for :class:`GKMMaturityMixin`."""


@deprecated("GKSSchemaMixin is deprecated; use GKMSchemaMixin instead.")
class GKSSchemaMixin(GKMSchemaMixin):
    """Deprecated alias for :class:`GKMSchemaMixin`."""


@deprecated("GKSMetadataMixin is deprecated; use GKMMetadataMixin instead.")
class GKSMetadataMixin(GKMMetadataMixin):
    """Deprecated alias for :class:`GKMMetadataMixin`."""
