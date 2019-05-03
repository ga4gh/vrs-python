"""vmc models, defined at runtime from the spec

**This file should not be imported directly.**

Users should `import vmc` and then `vmc.models` to access models.

This module reads the spec and generates classes at runtime.  The
advantage of this approach over models defined in code is that the
models are always in sync with the spec.

"""

import os

import pkg_resources
import python_jsonschema_objects as pjs

schema_dir = pkg_resources.resource_filename(__name__, "_data/schema")
schema_path = schema_dir + "/vmc.json"
schema_file = os.path.basename(schema_path)

classes = models = None

def _build_classes():
    """load/reload models"""
    global classes, models
    builder = pjs.ObjectBuilder(schema_path)
    classes = models = builder.build_classes(standardize_names=False)  # TODO: named_only=True, 
    return classes

_build_classes()
