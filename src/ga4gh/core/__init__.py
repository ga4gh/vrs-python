__all__ = """ga4gh_digest""".split()


import warnings
from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:    # pragma: nocover
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound

from ._internal.ga4gh_digest import ga4gh_digest

