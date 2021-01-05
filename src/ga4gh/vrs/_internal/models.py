"""Generate VRS models at runtime from the json schema

**This module should not be imported directly.**

Instead, users should use one of the following:

  * `from ga4gh.vrs import models`, and refer to models with the
    abbreviated name, e.g., `models.Allele` (recommended)

  * `import ga4gh.vrs`, and refer to models using the fully-qualified
    module name, e.g., `ga4gh.vrs.models.Allele`


This module reads the spec and generates classes at runtime.  The
advantage of this approach over models defined in code is that the
models are always in sync with the spec.

"""

import logging
import os

import pkg_resources

from ga4gh.core import build_models, build_class_referable_attribute_map


_logger = logging.getLogger(__name__)


if "VR_SCHEMA_DIR" in os.environ:
    _logger.warning("VR_SCHEMA_DIR is defined but being ignored; Use VRS_SCHEMA_DIR instead")

try:
    # specify VRS_SCHEMA_DIR to use a schema other than the one
    # embedded in vrs-python
    schema_dir = os.environ["VRS_SCHEMA_DIR"]
except KeyError:
    schema_dir = pkg_resources.resource_filename(__name__, "data/schema")

schema_path = schema_dir + "/vr.json"


models = None
class_refatt_map = None

def _load_vrs_models():
    """load/reload models from `schema_path`

    This function facilitates reloading changes to the schema during
    development.

    """

    global class_refatt_map, models
    models = build_models(schema_path, standardize_names=False)
    class_refatt_map = build_class_referable_attribute_map(models)
    return models


_load_vrs_models()
