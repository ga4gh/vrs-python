from typing import List, Literal, Optional
import json

from pydantic import Field, BaseModel, RootModel
from pydantic.json_schema import GenerateJsonSchema

from ga4gh.core import entity_models, domain_models
from ga4gh.vrs import models


MODULE_TO_REF = {
    "ga4gh.core.entity_models": "/ga4gh/schema/gks-common/1.x/core-im/json",
    "ga4gh.core.domain_models": "/ga4gh/schema/gks-common/1.x/domain-entities/json",
    "ga4gh.vrs.models": "/ga4gh/schema/vrs/2.x/json",
}


def create_model_module_map(*modules) -> dict[str, str]:
    """Creates a mapping from model names to their modules."""
    model_module_map = {}
    for module in modules:
        for attr_name in dir(module):
            model = getattr(module, attr_name)
            if (
                isinstance(model, type)
                and issubclass(model, (BaseModel, RootModel))
                and model.__module__ == module.__name__
            ):
                model_module_map[attr_name] = MODULE_TO_REF[model.__module__]
    return model_module_map


MODEL_REF_MAP = create_model_module_map(domain_models, entity_models, models)


class Allele(models.Allele, extra="forbid"):
    """The state of a molecule at a Location."""

    class Config:
        @staticmethod
        def json_schema_extra(cls):
            cls["properties"]["location"]["oneOf"] = cls["properties"]["location"][
                "anyOf"
            ]
            del cls["properties"]["location"]["anyOf"]

    expressions: Optional[List[models.Expression]] = Field(None, ordered=False)
    maturity: Literal["draft"] = Field("draft", frozen=True)
    alternativeLabels: Optional[List[str]] = Field(
        None, description="Alternative name(s) for the Entity.", ordered=False
    )
    extensions: Optional[List[entity_models.Extension]] = Field(
        None,
        description="A list of extensions to the Entity, that allow for capture of information not directly supported by elements defined in the model.",
        ordered=False,
    )


class GksGenerateJsonSchema(GenerateJsonSchema):
    def traverse_and_modify(self, schema):
        if isinstance(schema, dict):
            if "anyOf" in schema:
                schema["anyOf"] = [
                    s for s in schema["anyOf"] if s.get("type") != "null"
                ]
                if len(schema["anyOf"]) == 1:
                    schema.update(schema.pop("anyOf")[0])

            schema.pop("default", None)

            enum = schema.get("enum") or []
            if len(enum) == 1:
                del schema["enum"]

            if "properties" in schema:
                for prop in schema["properties"].values():
                    prop.pop("title", None)

            if "$ref" in schema:
                class_name = schema["$ref"].split("/")[-1]
                schema["$ref"] = f"{MODEL_REF_MAP[class_name]}/{class_name}"

            for value in schema.values():
                self.traverse_and_modify(value)

        elif isinstance(schema, list):
            for item in schema:
                self.traverse_and_modify(item)

    def generate(self, schema, mode="validation"):
        json_schema = super().generate(schema, mode=mode)
        json_schema["$schema"] = self.schema_dialect

        if "$defs" in json_schema:
            del json_schema["$defs"]

        model_class = schema.get("schema").get("cls")
        json_schema["title"] = model_class.__name__

        json_schema["properties"].pop("maturity", None)

        if "maturity" in model_class.model_fields:
            json_schema["maturity"] = model_class.model_fields["maturity"].default

        if "type" in model_class.model_fields:
            if "required" in json_schema:
                json_schema["required"].append("type")
            else:
                json_schema["required"] = ["type"]

        if hasattr(model_class, "is_ga4gh_identifiable"):
            if model_class.is_ga4gh_identifiable():
                json_schema["ga4ghDigest"] = {
                    "prefix": model_class.ga4gh.prefix,
                    "keys": model_class.ga4gh.keys,
                }

        json_schema["$id"] = (
            f"https://w3id.org/ga4gh/schema/vrs/2.x/json/{model_class.__name__}"
        )
        json_schema["description"] = model_class.__doc__

        self.traverse_and_modify(json_schema)
        return json_schema


with open("examples/Allele.json", "w") as wf:
    json.dump(
        Allele.model_json_schema(schema_generator=GksGenerateJsonSchema), wf, indent=2
    )
