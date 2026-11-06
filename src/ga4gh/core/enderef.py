"""Convert ga4gh objects between inlined and referenced forms.

Referable attributes are those that may be either an object or
CURIE.

These functions require a "class referable attribute" map, which will
typically be generated from the schema by
pydantic_class_refatt_map() in .models.py.

"""

import logging
from collections.abc import Mapping, MutableMapping, Sequence

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
    o: BaseModel,
    cra_map: Mapping[str, Sequence[str]],
    object_store: MutableMapping[str, BaseModel] | None = None,
    return_id_obj_tuple: bool = False,
) -> BaseModel | tuple[str | None, BaseModel]:
    """Recursively convert "referable attributes" from inlined to referenced form.  Returns a new object.

    :param o: VRS object
    :param cra_map: class referrable-attribute map; { o.type: [attr, attr, ...] }
    :param object_store: `collections.abc.MutableMapping` subclass that stores
        objects as they are referenced. If not given, referenced objects are not stored.
    :param return_id_obj_tuple: include ID with return value or not
    :return: new, enref'd object, as well as object ID if param set
    :raises TypeError: If a root-model CURIE reference has a non-GA4GH prefix.
        Plain string references are left unchanged.
    :raises ValueError: If the input is not an identifiable Pydantic model.
    """

    def _id_and_store(o: BaseModel) -> str | None:
        """Compute an identifier and store the object when a store is available.

        :param o: Object to process.
        :return: Computed identifier, or None for an unidentifiable object.
        """
        _id = ga4gh_identify(o)
        if _id and object_store is not None:
            object_store[_id] = o
        return _id

    def _enref(o: BaseModel) -> str | None:
        """Replace referable attributes in place and identify the visited object.

        :param o: Object to process.
        :return: Identifier of the updated object, or None.
        :raises TypeError: If a root-model CURIE has a non-GA4GH prefix.
        """
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


def ga4gh_deref(
    o: BaseModel,
    cra_map: Mapping[str, Sequence[str]],
    object_store: Mapping[str, BaseModel],
) -> BaseModel:
    """Copy an object and replace referable identifiers with inline objects.

    `object_store` must be a mappable object and is required for
    dereferencing.

    The top-level input is copied. Objects retrieved from the store are inserted
    directly and may be mutated while their own references are resolved.

    :param o: Object to process.
    :param cra_map: Mapping from model names to referable attribute names.
    :param object_store: Mapping from identifiers to referenced objects.
    :return: Copy with identifier references replaced by objects.
    :raises ValueError: If the input is not an identifiable Pydantic model.
    :raises KeyError: If a referenced object is absent from the store.
    """

    def _deref(o: BaseModel) -> BaseModel:
        """Resolve references recursively, mutating each visited object.

        :param o: Object to process.
        :return: Updated object with inline references.
        """
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
