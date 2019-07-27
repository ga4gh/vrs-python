import base64
import hashlib


def sha512t24u(blob):
    """generate a base64url-encode, truncated SHA-512 digest for given
    binary data

    The sha512t24u digest is a convention for constructing and
    formatting digests for use as object identifiers. Specifically::
    
        * generate a SHA512 digest on binary data
        * truncate at 24 bytes
        * encode using base64url encoding

    Examples:
    >>> sha512t24u(b'')
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    >>> sha512t24u(b"ACGT")
    'aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'

    """

    digest_size = 24
    digest = hashlib.sha512(blob).digest()
    tdigest_b64us = base64.urlsafe_b64encode(digest[:digest_size])
    return tdigest_b64us.decode("ascii")
