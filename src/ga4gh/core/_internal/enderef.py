"""Convert ga4gh objects between inlined and referenced forms.

Referable attributes are those that may be either an object or
CURIE.

These functions require a "class referable attribute" map, which will
typically be generated from the schema by
build_class_referable_attribute_map() in .models.py.

"""


from .identifiers import ga4gh_identify
from .jsonschema import is_identifiable, is_pjs_instance

def _ga4gh_copy(o):
    return o.__class__(**o.as_dict())


def ga4gh_enref(o, cra_map, object_store=None, max_depth=None):
    """Convert "referable attributes" in-place from inlined to referenced
    form.
    
    If `object_store` is provided, it must be a Mapping subclass and
    will be used to store objects as they are referenced.

    """

    assert is_pjs_instance(o)

    o = _ga4gh_copy(o)

    if o.type not in cra_map:
        return o
    
    for att in cra_map[o.type]:
        if not is_identifiable(o[att]):  # refatt is already a ref
            continue
        _id = ga4gh_identify(o[att])
        if object_store is not None:
            object_store[_id] = o[att]
        o[att] = _id
    return o


def ga4gh_deref(o, cra_map, object_store, max_depth=None):
    """Convert "referable attributes" in-place from referenced to inlined
    form.

    `object_store` must be a mappable object and is required for
    dereferencing.

    """

    assert is_pjs_instance(o)

    if o.type not in cra_map:
        return o
    
    o = _ga4gh_copy(o)

    for att in cra_map[o.type]:
        if is_identifiable(o[att]):  # refatt is already an instance
            continue
        _id = o[att]
        o[att] = object_store[_id]

    return o
