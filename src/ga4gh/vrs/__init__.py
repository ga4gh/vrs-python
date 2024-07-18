# flake8: noqa
"""Public interface to the GA4GH Variation Representation reference
implementation
"""
from importlib.metadata import version, PackageNotFoundError

from .normalize import normalize
from .enderef import vrs_deref, vrs_enref
from .models import VrsType
from . import models

__all__ = [
    "normalize",
    "vrs_deref",
    "vrs_enref",
    "VrsType",
    "models"
]

try:
    __version__ = version(__name__)
except PackageNotFoundError:    # pragma: nocover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
