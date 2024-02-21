"""Python support used across GA4GH projects

"""

import warnings
from importlib.metadata import version, PackageNotFoundError

from ._internal.digests import sha512t24u
from ._internal.enderef import ga4gh_enref, ga4gh_deref
from ._internal.exceptions import GA4GHError
from ._internal.identifiers import (
    ga4gh_digest, ga4gh_identify, ga4gh_serialize, is_ga4gh_identifier,
    parse_ga4gh_identifier, VrsObjectIdentifierIs, use_ga4gh_compute_identifier_when,
    CURIE_NAMESPACE, CURIE_SEP, GA4GH_PREFIX_SEP, GA4GH_IR_REGEXP, GA4GH_DIGEST_REGEXP
)
from ._internal.pydantic import (
    is_pydantic_instance, is_curie_type, is_ga4gh_identifiable, is_literal, pydantic_copy
)
from ._internal import models as core_models

try:
    __version__ = version(__name__)
except PackageNotFoundError:    # pragma: nocover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
