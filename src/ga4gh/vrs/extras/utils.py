"""miscellaneous utility functions for vrs.extras

"""

from base64 import urlsafe_b64decode, urlsafe_b64encode
from binascii import hexlify, unhexlify
import math


def _format_time(timespan, precision=3):
    """Formats the timespan in a human readable form

    >>> _format_time(0.35)
    '350 ms'

    >>> _format_time(35)
    '35 s'

    >>> _format_time(3500)
    '58min 20s'


    lovingly borrowed from
    https://github.com/ipython/ipython/blob/master/IPython/core/magics/execution.py

    """

    if timespan >= 60.0:
        # we have more than a minute, format that in a human readable form
        # Idea from http://snipplr.com/view/5713/
        parts = [("d", 60*60*24),("h", 60*60),("min", 60), ("s", 1)]
        time = []
        leftover = timespan
        for suffix, length in parts:
            value = int(leftover / length)
            if value > 0:
                leftover = leftover % length
                time.append(u"%s%s" % (str(value), suffix))
            if leftover < 1:
                break
        return " ".join(time)

    units = [u"s", u"ms", u"us", u"ns"]  # the save value
    scaling = [1, 1e3, 1e6, 1e9]

    if timespan > 0.0:
        order = min(-int(math.floor(math.log10(timespan)) // 3), 3)
    else:
        order = 3
    return u"%.*g %s" % (precision, timespan * scaling[order], units[order])


def hex_to_base64url(s):
    """convert hex string to base64 string"""
    return urlsafe_b64encode(unhexlify(s)).decode("ascii")

def base64url_to_hex(s):
    """convert base64 string to hex string"""
    return hexlify(urlsafe_b64decode(s)).decode("ascii")
