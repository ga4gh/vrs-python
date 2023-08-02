from ga4gh.core import ga4gh_deref, ga4gh_enref


def vrs_enref(o, object_store=None):
    raise RuntimeError("Not Implemented!")
    return ga4gh_enref(o, cra_map=class_refatt_map, object_store=object_store)


def vrs_deref(o, object_store):
    raise RuntimeError("Not Implemented!")
    return ga4gh_deref(o, cra_map=class_refatt_map, object_store=object_store)
