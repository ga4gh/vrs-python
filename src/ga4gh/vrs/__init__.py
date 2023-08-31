# flake8: noqa
"""Public interface to the GA4GH Variation Representation reference
implementation
"""
from importlib.metadata import version, PackageNotFoundError

from .normalize import normalize
from ._internal.enderef import vrs_deref, vrs_enref
from ._internal import models

__all__ = """models normalize schema_path vrs_deref vrs_enref""".split()

try:
    __version__ = version(__name__)
except PackageNotFoundError:    # pragma: nocover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
