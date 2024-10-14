import re


ref_re = re.compile(r":ref:`(.*?)(\s?<.*>)?`")
link_re = re.compile(r"`(.*?)\s?\<(.*)\>`_")


def scrub_rst_markup(string):
    string = ref_re.sub(r"\g<1>", string)
    string = link_re.sub(r"[\g<1>](\g<2>)", string)
    string = string.replace("\n", " ")
    return string
