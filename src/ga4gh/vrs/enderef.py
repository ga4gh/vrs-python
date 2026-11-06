"""Convert VRS objects between inline objects and identifier references.

Given an Allele with an inline SequenceLocation::

    from ga4gh.vrs import vrs_deref, vrs_enref

    store = {}
    allele_id, referenced = vrs_enref(allele, store, return_id_obj_tuple=True)
    restored = vrs_deref(referenced, store)

The round trip is::

    allele.location       -> SequenceLocation(...)
    referenced.location   -> "ga4gh:SL.<digest>"
    store[allele_id]       -> referenced Allele
    restored.location     -> stored SequenceLocation

Both operations copy the top-level input. Dereferencing inserts store values
directly, so resolving nested references can mutate those stored objects.
"""

from collections.abc import Mapping, MutableMapping

from pydantic.main import BaseModel

from ga4gh.core import ga4gh_deref, ga4gh_enref

from .models import class_refatt_map


def vrs_enref(
    o: BaseModel,
    object_store: MutableMapping[str, BaseModel] | None = None,
    return_id_obj_tuple: bool = False,
) -> BaseModel | tuple[str | None, BaseModel]:
    """Copy a VRS object and replace referable attributes with identifiers.

    :param o: Identifiable VRS object to copy.
    :param object_store: Optional store for objects replaced by references.
    :param return_id_obj_tuple: Include the identifier alongside the copied object.
    :return: Referenced copy, optionally paired with its identifier.
    :raises ValueError: If the input is not an identifiable Pydantic object.
    :raises TypeError: If a root-model CURIE reference has a non-GA4GH prefix.
        Plain string references are left unchanged.
    """
    return ga4gh_enref(
        o,
        cra_map=class_refatt_map,
        object_store=object_store,
        return_id_obj_tuple=return_id_obj_tuple,
    )


def vrs_deref(o: BaseModel, object_store: Mapping[str, BaseModel]) -> BaseModel:
    """Copy a VRS object and resolve its identifier references.

    :param o: Identifiable VRS object to copy.
    :param object_store: Store mapping identifiers to referenced objects. Retrieved
        objects are inserted directly and may be mutated during dereferencing.
    :return: Copy with references replaced by inline objects.
    :raises ValueError: If the input is not an identifiable Pydantic object.
    :raises KeyError: If a referenced identifier is absent from the store.
    """
    return ga4gh_deref(o, cra_map=class_refatt_map, object_store=object_store)
