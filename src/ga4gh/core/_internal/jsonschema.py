"""Support for JSON Schema and GA4GH schema conventions using
python-jsonschema-objects

This module has evolved to include several different kinds of
functions:

 * schema functions operate on entire json schemas and
   subschemas
 * pjs class functions interrogate python_jsonschema_objects (pjs) classes
 * instance functions interrogate pjs objects/instances
 * attribute functions interrogate pjs attributes

"""

import logging

import python_jsonschema_objects as pjs
import yaml

from .helpers import pjs_filter

_logger = logging.getLogger(__name__)


############################################################################
# Schema Functions

def build_models(path, standardize_names=False):
    """load models from json schema at path"""
    with open(path, "r") as yaml_file:
        y = yaml.load(yaml_file, yaml.SafeLoader)
    builder = pjs.ObjectBuilder(pjs_filter(y))
    models = builder.build_classes(standardize_names=standardize_names)
    return models


def build_class_referable_attribute_map(models):
    """given a set of json schema models (as from build_models()), return
    a map of classes ⇒ referrable attributes

    The practical meaning is to define a map like
       {"Allele": ["location"], ...}
    for any json schema class that has one or more attributes
    that may be either an object or a CURIE.

    """
    cra_map = {cn: get_referable_attributes(models[cn]) for cn in models}
    return {cn: refatts for cn, refatts in cra_map.items()
            if refatts}


def get_referable_attributes(cls):
    """for a given pjs class, return list of attributes that may be either
    objects or CURIEs

    ie, a list of attributes that may be inlined objects or references
    to them

    """

    if not is_pjs_class(cls):
        return None
    atts = cls.__prop_names__
    refatts = [att for att in atts if is_referable(cls.propinfo(att))]
    return refatts


def is_referable(json_subschema):
    """return True if field is a referable object, or a list of referable
    objects

    """

    if "oneOf" in json_subschema:
        # schema is oneOf of a CURIE and non-CURIE type
        refs = [oo.get("$ref", None) for oo in json_subschema["oneOf"]]
        return (any(r and r.endswith("/CURIE") for r in refs)
                    and any(not r.endswith("/CURIE") for r in refs))

    if "type" in json_subschema:
        t = json_subschema["type"]
        if t == "array":
            # an array of referable types
            return is_referable(json_subschema["items"])

    return False


############################################################################
# Class Functions
# (argument is a pjs class)

def is_pjs_class(c):
    """return True if argument is a pjs class object

    """
    mro = getattr(c, "__mro__", [])
    return pjs.classbuilder.ProtocolBase in mro


############################################################################
# Instance/object Functions
# (argument is a pjs object instance)

def is_pjs_instance(o):
    """return True if object is a python jsonschema object"""
    return isinstance(o, pjs.classbuilder.ProtocolBase)


def is_pjs_literal(o):
    """return True if object is a python jsonschema object literal"""
    return isinstance(o, pjs.literals.LiteralValue)


def is_pjs_array(o):
    """return True if object is a python jsonschema object array"""
    return getattr(o, "type", None) == "array"


def is_curie_type(o):
    """return True if object is a python jsonschema class that represents
    a CURIE, e.g., sequence_id

    """
    return o.__class__.__name__.endswith("/CURIE")


def pjs_copy(o):
    """create a new instance of a pjs object.

    A bug in pjs prevents using Python's more efficient copy module.

    """
    return o.__class__(**o.as_dict())


# backward compatibility
is_literal = is_pjs_literal
is_array = is_pjs_array
is_curie = is_curie_type


############################################################################
# Attribute Functions

def is_identifiable(o):
    """return True if object is identifiable

    An object is considered identifiable if it contains an `_id`
    attribute

    """

    return is_pjs_instance(o) and ("_id" in o)
