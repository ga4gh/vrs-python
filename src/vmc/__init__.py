# flake8: noqa
"""vmc -- public interface to vmcdemo

Most callers should use something like this:

>>> from vmc import models, computed_id, digest, serialize

"""

from ._models import models, schema_path
from .digest import computed_id, computed_identifier, digest, serialize
from .serialize import serialize
#from .utils import get_vmc_sequence_identifier, ir_to_id, id_to_ir
