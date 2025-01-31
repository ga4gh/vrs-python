"""Python support used across GA4GH projects"""

from importlib.metadata import version, PackageNotFoundError

from ga4gh.core.digests import sha512t24u
from ga4gh.core.enderef import ga4gh_enref, ga4gh_deref
from ga4gh.core.identifiers import (
    ga4gh_digest,
    ga4gh_identify,
    ga4gh_serialize,
    is_ga4gh_identifier,
    VrsObjectIdentifierIs,
    use_ga4gh_compute_identifier_when,
    CURIE_NAMESPACE,
    CURIE_SEP,
    GA4GH_PREFIX_SEP,
    GA4GH_IR_REGEXP,
    GA4GH_DIGEST_REGEXP,
    PrevVrsVersion,
)
from ga4gh.core.pydantic import is_pydantic_instance, is_curie_type, pydantic_copy
from ga4gh.core import models as core_models

__all__ = [
    "CURIE_NAMESPACE",
    "CURIE_SEP",
    "GA4GH_DIGEST_REGEXP",
    "GA4GH_IR_REGEXP",
    "GA4GH_PREFIX_SEP",
    "PrevVrsVersion",
    "VrsObjectIdentifierIs",
    "core_models",
    "ga4gh_deref",
    "ga4gh_digest",
    "ga4gh_enref",
    "ga4gh_identify",
    "ga4gh_serialize",
    "is_curie_type",
    "is_ga4gh_identifier",
    "is_pydantic_instance",
    "pydantic_copy",
    "sha512t24u",
    "use_ga4gh_compute_identifier_when",
]

try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: nocover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
