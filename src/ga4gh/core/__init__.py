"""Python support used across GA4GH projects

"""

import warnings
from importlib.metadata import version, PackageNotFoundError

from ._internal.digests import sha512t24u
from ._internal.enderef import ga4gh_enref, ga4gh_deref
from ._internal.exceptions import GA4GHError
from ._internal.identifiers import (
    ga4gh_digest, ga4gh_identify, ga4gh_serialize, is_ga4gh_identifier,
    parse_ga4gh_identifier
)
from ._internal.pydantic import (
    is_pydantic_instance, is_curie_type, is_identifiable, is_literal, pydantic_copy
)

try:
    __version__ = version(__name__)
except PackageNotFoundError:    # pragma: nocover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
