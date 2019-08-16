"""Support for JSON Schema and GA4GH schema conventions 

"""

import python_jsonschema_objects as pjs


def build_models(path, standardize_names):
    """load models from json schema at path"""
    builder = pjs.ObjectBuilder(path)
    models = builder.build_classes(standardize_names=False)
    return models


def ga4gh_as_dict(o, filter_optional_keys=True):
    """return the VR objects as a dictionary, optionally filtering optional keys"""

    def filter_dict(d):
        try:
            return {k: filter_dict(d[k])
                    for k in d
                    if not k.startswith("_")}
        except:
            return d


    if not is_class(o):
        raise ValueError("Attempted to call ga4gh_as_dict on non-VR object")
    
    d = o.as_dict()

    if filter_optional_keys:
        d = filter_dict(d)

    return d


def is_class(o):
    """return True if object is a python jsonschema object"""
    return isinstance(o, pjs.classbuilder.ProtocolBase)


def is_identifiable(o):
    """return True if object is identifiable

    An object is considered identifiable if it contains an `_id` attribute
    """
    return is_class(o) and ("_id" in o)


def is_literal(o):
    """return True if object is a python jsonschema object literal"""
    return isinstance(o, pjs.literals.LiteralValue)


