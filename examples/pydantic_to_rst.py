from enum import Enum
from annotated_types import MaxLen, MinLen
from pydantic import BaseModel
from pydantic._internal._model_construction import ModelMetaclass
from typing import Annotated, Any, List
from typing import get_args, get_origin

from utils import EXCLUDE_PROPS
from pydantic_to_json_schema import Allele


PYTHON_TO_JSON_TYPES = {
    str: "string",
    int: "integer",
    float: "float",
    None: "null",
    list: "array",
    bool: "boolean",
    dict: "object",
}


def get_field_type(field_anno: Any) -> str | None:
    if field_anno in PYTHON_TO_JSON_TYPES:
        field_type = PYTHON_TO_JSON_TYPES[field_anno]
    elif isinstance(field_anno, ModelMetaclass):
        field_type = f":ref:`{field_anno.__name__}`"
    elif get_origin(field_anno) is Annotated:
        field_type = "string"
    elif isinstance(field_anno, Enum):
        field_type = field_anno.name
    else:
        field_type = None
    return field_type


def get_limits(
    field_name: str, field_metadata: list[MinLen | MaxLen], field_is_list: bool
) -> list[str]:
    min_length = None
    max_length = None

    if field_metadata:
        if len(field_metadata) == 2:
            min_length = field_metadata[0].min_length
            max_length = field_metadata[1].max_length
        else:
            min_length = getattr(field_metadata[0], "min_length", None)
            max_length = getattr(field_metadata[0], "max_length", None)

    if min_length is not None:
        limits = [min_length]
    else:
        limits = [1 if field_name == "type" else 0]

    if max_length is not None:
        limits.append(max_length)
    else:
        limits.append("m" if field_is_list else 1)
    return [str(limit) for limit in limits]


def generate(model: BaseModel) -> str:
    # TODO: Inheritance
    inheritance = ""

    rst_data = [
        "**Computational Definition**",
        "",
        model.__doc__,
        "",
        "**Information Model**",
        "",
        inheritance,
        ".. list-table::",
        "   :class: clean-wrap",
        "   :header-rows: 1",
        "   :align: left",
        "   :widths: auto",
        "",
        "   * - Field",
        "     - Type",
        "     - Limits",
        "     - Description",
    ]

    for field_name, field_info in model.model_fields.items():
        if field_name in EXCLUDE_PROPS:
            continue

        field_type = None
        if isinstance(field_info.annotation, type) and issubclass(
            field_info.annotation, Enum
        ):
            field_is_list = False
            field_type = f":ref:`{field_info.annotation.__name__}`"
        else:
            field_annotation = tuple(
                anno
                for anno in get_args(field_info.annotation)
                if anno is not type(None)
            )

            field_is_list = get_origin(field_annotation[0]) in {list, List}

        limits = get_limits(field_name, field_info.metadata, field_is_list)

        if field_type is None:
            if len(field_annotation) > 1:
                field_type = " | ".join(
                    [get_field_type(anno) for anno in field_annotation]
                )
            else:
                field_type = get_field_type(
                    get_args(field_annotation[0])[0]
                    if field_is_list
                    else field_annotation[0]
                )

        rst_data.extend(
            [
                f"   * - {field_name}",
                f"     - {field_type}",
                f"     - {'..'.join(limits)}",
                f"     - {field_info.description}",
            ]
        )

    return "\n".join(rst_data)


if __name__ == "__main__":
    with open("examples/Allele.rst", "w") as wf:
        wf.write(generate(Allele))
