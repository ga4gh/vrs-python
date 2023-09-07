from ga4gh.core import ga4gh_deref, ga4gh_enref

from .models import class_refatt_map


def vrs_enref(o, object_store=None):
    return ga4gh_enref(o, cra_map=class_refatt_map, object_store=object_store)


def vrs_deref(o, object_store):
    return ga4gh_deref(o, cra_map=class_refatt_map, object_store=object_store)
