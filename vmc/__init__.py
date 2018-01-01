# flake8: noqa
"""vmc -- public interface to vmcdemo

Most callers should use something like this:

>>> from vmc import models, computed_id, vmc_digest, vmc_serialize

"""

from ._models import models, schema_path
from .digest import computed_id, computed_identifier, vmc_digest, vmc_serialize
from .seqrepo import get_vmc_sequence_id
