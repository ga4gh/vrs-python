# flake8: noqa
"""Public interface to the GA4GH Variation Representation reference
implementation

"""

from pkg_resources import get_distribution, DistributionNotFound

from .normalize import normalize
from ._internal.enderef import vrs_deref, vrs_enref
from ._internal import models
schema_path = models.schema_path

__all__ = """models normalize schema_path vrs_deref vrs_enref""".split()

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:    # pragma: nocover
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound
