"""Public interface to the GA4GH Variation Representation reference implementation"""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as package_version

from ga4gh.vrs import models
from ga4gh.vrs.enderef import vrs_deref, vrs_enref
from ga4gh.vrs.models import VrsType
from ga4gh.vrs.normalize import normalize
from ga4gh.vrs.version import VRS_VERSION

try:
    __version__ = package_version(__name__)
except PackageNotFoundError:  # pragma: nocover
    __version__ = "unknown"
finally:
    del package_version, PackageNotFoundError


__all__ = [
    "VRS_VERSION",
    "VrsType",
    "models",
    "normalize",
    "vrs_deref",
    "vrs_enref",
]
