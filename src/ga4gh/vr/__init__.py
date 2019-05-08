# flake8: noqa
"""Public interface to the GA4GH Variation Representation reference
implementation

Import as:

>>> from ga4gh.vr import models, computed_id, digest, serialize

Deprecation and change notices are provided only for definitions
obtained by importing ga4gh.vr as shown above.

Modules names that begin with an underscore are internal; definitions
therein are not part of the public interface and may change or
disappear without notice.

"""

import warnings
from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound


from ._models import models, schema_path
from .digest import computed_id, computed_identifier, digest, serialize
from .serialize import serialize
#from .utils import get_vmc_sequence_identifier, ir_to_id, id_to_ir

