"""Python support used across GA4GH projects"""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as package_version

import ga4gh.core.models as core_models
from ga4gh.core.digests import sha512t24u
from ga4gh.core.enderef import ga4gh_deref, ga4gh_enref
from ga4gh.core.identifiers import (
    CURIE_NAMESPACE,
    CURIE_SEP,
    GA4GH_DIGEST_REGEXP,
    GA4GH_IR_REGEXP,
    GA4GH_PREFIX_SEP,
    PrevVrsVersion,
    VrsObjectIdentifierIs,
    ga4gh_digest,
    ga4gh_identify,
    ga4gh_serialize,
    is_ga4gh_identifier,
    use_ga4gh_compute_identifier_when,
)
from ga4gh.core.metadata import (
    GKSMaturityMixin,
    GKSMetadataMixin,
    GKSSchemaMixin,
    Maturity,
)
from ga4gh.core.models import GKSCoreMetadataMixin
from ga4gh.core.pydantic import is_curie_type, is_pydantic_instance, pydantic_copy
from ga4gh.core.version import CORE_VERSION

try:
    __version__ = package_version(__name__)
except PackageNotFoundError:  # pragma: nocover
    __version__ = "unknown"
finally:
    del package_version, PackageNotFoundError


__all__ = [
    "CORE_VERSION",
    "CURIE_NAMESPACE",
    "CURIE_SEP",
    "GA4GH_DIGEST_REGEXP",
    "GA4GH_IR_REGEXP",
    "GA4GH_PREFIX_SEP",
    "GKSCoreMetadataMixin",
    "GKSMaturityMixin",
    "GKSMetadataMixin",
    "GKSSchemaMixin",
    "Maturity",
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
