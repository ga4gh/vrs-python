"""Generate VR models at runtime from the spec

**This file should not be imported directly.**

Users should use one of the following:

  * `from ga4gh.vr import models`, and refer to models with the
    abbreviated name, e.g., `models.Allele` (recommended)

  * `import ga4gh.vr`, and refer to models using the fully-qualified
    module name, e.g., `ga4gh.vr.models.Allele`


This module reads the spec and generates classes at runtime.  The
advantage of this approach over models defined in code is that the
models are always in sync with the spec.

"""

import os

import pkg_resources

from ga4gh.core import build_classes


schema_dir = pkg_resources.resource_filename(__name__, "data/schema")
schema_path = schema_dir + "/vr.json"
schema_file = os.path.basename(schema_path)

classes = models = None

def _build_classes():
    """load/reload models
    
    developers may call this function to reload schemas during development
    """
    global classes, models
    classes = models = build_classes(schema_path, standardize_names=False)
    return classes

_build_classes()
