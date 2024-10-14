import re


REF_RE = re.compile(r":ref:`(.*?)(\s?<.*>)?`")
LINK_RE = re.compile(r"`(.*?)\s?\<(.*)\>`_")
EXCLUDE_PROPS = {"maturity"}


def scrub_rst_markup(string):
    string = REF_RE.sub(r"\g<1>", string)
    string = LINK_RE.sub(r"[\g<1>](\g<2>)", string)
    string = string.replace("\n", " ")
    return string
