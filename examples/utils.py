import re
from typing import get_args

from pydantic import RootModel

from ga4gh.core import entity_models, domain_models
from ga4gh.vrs import models


REF_RE = re.compile(r":ref:`(.*?)(\s?<.*>)?`")
LINK_RE = re.compile(r"`(.*?)\s?\<(.*)\>`_")
EXCLUDE_PROPS = {"maturity"}


def scrub_rst_markup(string):
    string = REF_RE.sub(r"\g<1>", string)
    string = LINK_RE.sub(r"[\g<1>](\g<2>)", string)
    string = string.replace("\n", " ")
    return string


def map_abc_to_instances(*modules) -> dict[str, str]:
    """Creates a mapping from ABC model names to their instance model names."""
    abc_to_instances = {}
    excluded_types = {
        "list",
        "str",
        "Optional",
        "int",
        "float",
        "dict",
        "bool",
        "set",
        "tuple",
    }

    for module in modules:
        for attr_name in dir(module):
            model = getattr(module, attr_name)
            if (
                isinstance(model, type)
                and issubclass(model, RootModel)
                and model.__module__ == module.__name__
            ):
                root_anno = model.model_fields["root"].annotation
                root_annos = get_args(root_anno) or (root_anno,)
                root_anno_cls_names = [
                    cls.__name__
                    for cls in root_annos
                    if cls.__name__ not in excluded_types
                ]
                if root_anno_cls_names:
                    abc_to_instances[model.__name__] = root_anno_cls_names

    return abc_to_instances


def get_abc(key) -> str:
    """Get original ABC class name

    :param key: Class name
    :return: Original ABC class name
    """
    while key in INSTANCE_TO_ABC:
        key = INSTANCE_TO_ABC[key]
    return key


ABC_TO_INSTANCES = map_abc_to_instances(models, entity_models, domain_models)
INSTANCE_TO_ABC = {
    value: key for key, value_list in ABC_TO_INSTANCES.items() for value in value_list
}
INSTANCE_TO_ABC = {key: get_abc(key) for key in INSTANCE_TO_ABC}
