# flake8: noqa
"""Public interface to the GA4GH Variation Representation reference
implementation

"""


__all__ = """models normalize schema_path vr_deref vr_enref""".split()



import warnings
from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:    # pragma: nocover
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound


from .normalize import normalize
from ._internal.enderef import vr_deref, vr_enref
from ._internal.models import models, schema_path
