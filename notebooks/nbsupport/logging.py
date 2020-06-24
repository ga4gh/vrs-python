import logging

def set_log_level(level="WARNING"):
    logger = logging.getLogger()
    logger.setLevel(level)
    for h in logger.handlers:
        h.setLevel(level)
