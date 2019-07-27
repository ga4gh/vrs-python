__all__ = ("sha512t24u", "GA4GHError")


import warnings
from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:    # pragma: nocover
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound

from ._internal.exceptions import GA4GHError
from ._internal.digests import sha512t24u
