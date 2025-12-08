"""Convert ga4gh objects between inlined and referenced forms.

Referable attributes are those that may be either an object or
CURIE.

These functions require a "class referable attribute" map, which will
typically be generated from the schema by
build_class_referable_attribute_map() in .models.py.

"""

import logging

from pydantic.main import BaseModel

from .identifiers import ga4gh_identify, is_ga4gh_identifier
from .pydantic import (
    get_pydantic_root,
    is_curie_type,
    is_pydantic_instance,
    pydantic_copy,
)

_logger = logging.getLogger(__name__)


def ga4gh_enref(
    o,  # noqa: ANN001
    cra_map,  # noqa: ANN001
    object_store=None,  # noqa: ANN001
    return_id_obj_tuple: bool = False,
) -> tuple:
    """Recursively convert "referable attributes" from inlined to
    referenced form.  Returns a new object.

    :param o: VRS object
    :param cra_map: class referrable-attribute map; { o.type: [attr, attr, ...] }
    :param object_store: `collections.abc.MutableMapping` subclass that stores
        objects as they are referenced. If not given, referenced objects are not stored.
    :param return_id_obj_tuple: include ID with return value or not
    :return: new, enref'd object, as well as object ID if param set
    :raise TypeError: if any object IDs are non-GA4GH CURIEs
    """

    def _id_and_store(o) -> str | None:  # noqa: ANN001
        _id = ga4gh_identify(o)
        if _id and object_store is not None:
            object_store[_id] = o
        return _id

    def _enref(o: BaseModel) -> str | None:
        """depth-first recursive, in-place enref of object; returns id of object"""
        ref_att_names = cra_map.get(o.type, [])
        for ran in ref_att_names:
            v = getattr(o, ran)
            if isinstance(v, list):
                setattr(o, ran, [_enref(o2) for o2 in v])
            elif isinstance(v, str):
                pass
            elif is_curie_type(v):  # already a reference
                if not is_ga4gh_identifier(v):
                    msg = (
                        "Identifiable attribute CURIE is contains an invalid identifier"
                    )
                    raise TypeError(msg)
            elif v is not None:
                _id = _id_and_store(v)
                if _id:
                    setattr(o, ran, _id)

        return _id_and_store(o)

    if not is_pydantic_instance(o):
        msg = "Called ga4gh_enref() with non-pydantic instance"
        raise ValueError(msg)
    if not o.is_ga4gh_identifiable():
        msg = "Called ga4gh_enref() with non-identifiable object"
        raise ValueError(msg)

    # in-place replacement on object copy
    o = pydantic_copy(o)
    _id = _enref(o)
    return (_id, o) if return_id_obj_tuple else o


def ga4gh_deref(o, cra_map, object_store) -> BaseModel:  # noqa: ANN001
    """Convert "referable attributes" in-place from referenced to inlined
    form.

    `object_store` must be a mappable object and is required for
    dereferencing.

    Raises KeyError if any object cannot be dereferenced

    """

    def _deref(o: BaseModel):  # noqa: ANN202
        """depth-first recursive, in-place deref of object; returns id of object"""
        if o.type not in cra_map:
            _logger.warning("%s not in cra_map %s", o.type, cra_map)
            return o

        ref_att_names = cra_map[o.type]
        for ran in ref_att_names:
            v = getattr(o, ran)
            if isinstance(v, list):
                setattr(o, ran, [_deref(object_store[str(curie)]) for curie in v])
            elif is_ga4gh_identifier(v):
                v = get_pydantic_root(v)
                dereffed_identifier = object_store[str(v)]
                setattr(o, ran, _deref(dereffed_identifier))
            else:
                pass  # some object; pass as-is

        return o

    if not is_pydantic_instance(o):
        msg = "Called ga4gh_deref() with non-pydantic instance"
        raise ValueError(msg)
    if not o.is_ga4gh_identifiable():
        msg = "Called ga4gh_deref() with non-identifiable object"
        raise ValueError(msg)

    # in-place replacement on object copy
    o = pydantic_copy(o)
    _deref(o)
    return o
