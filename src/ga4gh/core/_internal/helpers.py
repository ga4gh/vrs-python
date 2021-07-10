"""miscellaneous functions for vrs-python

"""

import logging

_logger = logging.getLogger(__name__)


def pjs_filter(yaml_dict):
    """filter out schema elements that are not supported by python
    jsonschema objects (yet)

    """

    for message_name, message_definition in yaml_dict["definitions"].items():
        if "anyOf" in message_definition:
            _logger.warning("Removing anyOf attribute from %s to avoid python-jsonschema-objects error.", message_name)
            del message_definition["anyOf"]
        if "allOf" in message_definition:
            _logger.warning("Removing allOf attribute from %s to avoid python-jsonschema-objects error.", message_name)
            del message_definition["allOf"]
    return yaml_dict
