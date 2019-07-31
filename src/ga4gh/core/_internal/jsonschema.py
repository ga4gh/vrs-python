"""Provide support for JSON Schema using python-jsonschema-objects"""

import python_jsonschema_objects as pjs


def build_classes(path, standardize_names):
    """load/reload models"""
    builder = pjs.ObjectBuilder(path)
    classes = models = builder.build_classes(standardize_names=False)  # TODO: named_only=True, 
    return classes

def is_class(vro):
    """return True if object is a python jsonschema object"""
    return isinstance(vro, pjs.classbuilder.ProtocolBase)

def is_identifiable(vro):
    """return True if object is identifiable"""
    return is_class(vro) and ("_id" in vro)

def is_literal(vro):
    return isinstance(vro, pjs.literals.LiteralValue)
