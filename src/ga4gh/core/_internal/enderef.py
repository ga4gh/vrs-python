"""Convert ga4gh objects between inlined and referenced forms.

Referable attributes are those that may be either an object or
CURIE.

These functions require a "class referable attribute" map, which will
typically be generated from the schema by
build_class_referable_attribute_map() in .models.py.

"""

import logging

from .identifiers import ga4gh_identify, is_ga4gh_identifier
from .jsonschema import is_array, is_curie, is_identifiable, is_pjs_instance, pjs_copy

_logger = logging.getLogger(__name__)


def ga4gh_enref(o, cra_map, object_store=None):
    """Recursively convert "referable attributes" from inlined to
    referenced form.  Returns a new object.

    If `object_store` is provided, it must be an
    collections.abc.MutableMapping subclass and will be used to store
    objects as they are referenced.  If `object_store` is not provided
    (or None), referenced objects will not be stored.

    """
    def _id_and_store(o):
        _id = ga4gh_identify(o)
        if object_store is not None:
            object_store[_id] = o
        return _id

    def _enref(o):
        """depth-first recursive, in-place enref of object; returns id of object"""
        ref_att_names = cra_map.get(o.type, [])
        for ran in ref_att_names:
            v = o[ran]
            if is_array(v):
                o[ran] = [_enref(o2) for o2 in v]
            elif is_curie(v):    # already a reference
                assert is_ga4gh_identifier(v), "Identifiable attribute CURIE is contains an invalid identifier"
            elif v is not None:
                o[ran] = _id_and_store(v)

        return _id_and_store(o)

    if not is_pjs_instance(o):
        raise ValueError("Called ga4gh_enref() with non-python_jsonschema_object instance")
    if not is_identifiable(o):
        raise ValueError("Called ga4gh_enref() with non-identifiable object")

    # in-place replacement on object copy
    o = pjs_copy(o)
    _enref(o)
    return o


def ga4gh_deref(o, cra_map, object_store):
    """Convert "referable attributes" in-place from referenced to inlined
    form.

    `object_store` must be a mappable object and is required for
    dereferencing.

    Raises KeyError if any object cannot be dereferenced

    """
    def _deref(o):
        """depth-first recursive, in-place deref of object; returns id of object"""
        if o.type not in cra_map:
            return o

        ref_att_names = cra_map[o.type]
        for ran in ref_att_names:
            v = o[ran]
            if is_array(v):
                o[ran] = [_deref(object_store[str(curie)]) for curie in v]
            elif is_ga4gh_identifier(v):
                o[ran] = _deref(object_store[str(v)])
            else:
                pass    # some object; pass as-is

        return o

    if not is_pjs_instance(o):
        raise ValueError("Called ga4gh_deref() with non-python_jsonschema_object instance")
    if not is_identifiable(o):
        raise ValueError("Called ga4gh_deref() with non-identifiable object")

    # in-place replacement on object copy
    o = pjs_copy(o)
    _deref(o)
    return o
