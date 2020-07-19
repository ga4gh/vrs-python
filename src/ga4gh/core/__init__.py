import warnings
from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:    # pragma: nocover
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound


from ._internal.digests import sha512t24u
from ._internal.enderef import ga4gh_enref, ga4gh_deref
from ._internal.exceptions import GA4GHError
from ._internal.identifiers import ga4gh_digest, ga4gh_identify, ga4gh_serialize, is_ga4gh_identifier, parse_ga4gh_identifier
from ._internal.jsonschema import build_models, build_class_referable_attribute_map, is_pjs_class, is_pjs_instance, is_curie, is_identifiable, is_literal, pjs_copy
