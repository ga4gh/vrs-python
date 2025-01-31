"""Public interface to the GA4GH Variation Representation reference
implementation
"""

from importlib.metadata import PackageNotFoundError, version

from ga4gh.vrs import models
from ga4gh.vrs.enderef import vrs_deref, vrs_enref
from ga4gh.vrs.models import VrsType
from ga4gh.vrs.normalize import normalize

__all__ = ["VrsType", "models", "normalize", "vrs_deref", "vrs_enref"]

try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: nocover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
