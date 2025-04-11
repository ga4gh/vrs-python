"""Public interface to the GA4GH Variation Representation reference implementation"""

from importlib.metadata import PackageNotFoundError, version

from ga4gh.vrs import models
from ga4gh.vrs.enderef import vrs_deref, vrs_enref
from ga4gh.vrs.models import VrsType
from ga4gh.vrs.normalize import normalize

try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: nocover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError


VRS_VERSION = "2.0.1"

__all__ = [
    "VRS_VERSION",
    "VrsType",
    "models",
    "normalize",
    "vrs_deref",
    "vrs_enref",
]
