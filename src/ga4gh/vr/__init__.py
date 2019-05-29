# flake8: noqa
"""Public interface to the GA4GH Variation Representation reference
implementation

Import as:

>>> from ga4gh.vr import models, computed_id, serialize

Deprecation and change notices are provided only for definitions
obtained by importing ga4gh.vr as shown above.

Modules names that begin with an underscore are internal; definitions
therein are not part of the public interface and may change or
disappear without notice.

"""


__all__ = """models schema_path   compute_id ga4gh_digest serialize""".split()


import warnings
from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound


from ._internal.models import models, schema_path
from ._internal.ids import compute_id, ga4gh_digest, serialize
