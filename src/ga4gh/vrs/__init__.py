# flake8: noqa
"""Public interface to the GA4GH Variation Representation reference
implementation

"""


__all__ = """models normalize schema_path vrs_deref vrs_enref""".split()


from pkg_resources import get_distribution, DistributionNotFound

from ._internal.enderef import vrs_deref, vrs_enref
from ._internal.models import models, schema_path
from .normalize import normalize

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:    # pragma: nocover
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound
