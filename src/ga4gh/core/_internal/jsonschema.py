"""Support for JSON Schema and GA4GH schema conventions 

"""

import python_jsonschema_objects as pjs


def build_models(path, standardize_names):
    """load models from json schema at path"""
    builder = pjs.ObjectBuilder(path)
    models = builder.build_classes(standardize_names=False)
    return models


def is_class(o):
    """return True if object is a python jsonschema object"""
    return isinstance(o, pjs.classbuilder.ProtocolBase)


def is_curie(o):
    """return True if object is a python jsonschema class that represents
    a CURIE, e.g., sequence_id"""
    return o.__class__.__name__.endswith("/CURIE")


def is_identifiable(o):
    """return True if object is identifiable

    An object is considered identifiable if it contains an `_id` attribute
    """
    return is_class(o) and ("_id" in o)


def is_literal(o):
    """return True if object is a python jsonschema object literal"""
    return isinstance(o, pjs.literals.LiteralValue)


def is_array(o):
    """return True if object is a python jsonschema object array"""
    return getattr(o, "type", None) == "array"
