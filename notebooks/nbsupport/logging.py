import logging


def set_log_level(level):
    """reset log level for the root logger and all handlers"""
    logger = logging.getLogger()
    logger.setLevel(level)
    for h in logger.handlers:
        h.setLevel(level)
