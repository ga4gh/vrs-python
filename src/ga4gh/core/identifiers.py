"""serializes, digests, and identifies GA4GH objects

In GA4GH schemas with nested objects, serialize, digest, and identify
are entangled.

For example, here is a call path for ga4gh_identify called on an Allele:
    ga4gh_identify(allele)
    + ga4gh_digest(allele)
    ++ ga4gh_serialize(allele)
    +++ ga4gh_digest(allele.location)
    ++++ ga4gh_serialize(allele.location)
    +++ ga4gh_serialize(allele.state)

For that reason, they are implemented here in one file.

"""
from canonicaljson import encode_canonical_json
import contextvars
import re
from contextlib import ContextDecorator
from enum import Enum, IntEnum
from typing import Optional
from pydantic import BaseModel

from ga4gh.core.pydantic import get_pydantic_root

__all__ = "ga4gh_digest ga4gh_identify ga4gh_serialize is_ga4gh_identifier".split()

CURIE_NAMESPACE = "ga4gh"
CURIE_SEP = ":"
GA4GH_PREFIX_SEP = "."

GA4GH_IR_REGEXP = re.compile(r"^ga4gh:(?P<type>[^.]+)\.(?P<digest>[0-9A-Za-z_\-]{32})$")
GA4GH_DIGEST_REGEXP = re.compile(r"^[0-9A-Za-z_\-]{32}$")

NS_W_SEP = f"{CURIE_NAMESPACE}{CURIE_SEP}"


class VrsObjectIdentifierIs(IntEnum):
    """
    Defines the state for when the `ga4gh_identify` method should compute
    an identifier ('id' attribute) for the specified object.  The options are:
      ANY - Always compute the identifier (this is the default behavior)
      GA4GH_INVALID - Compute the identifier if it is missing or is present but syntactically invalid
      MISSING - Only compute the identifier if missing

    The default behavior is safe and ensures that the identifiers are correct,
    but at a performance cost. Where the source of inputs to `ga4gh_identify`
    are well controlled, for example when annotating a VCF file with VRS IDs,
    using `MISSING` can improve performance.
    """

    ANY = 0
    GA4GH_INVALID = 1
    MISSING = 2


class PrevVrsVersion(str, Enum):
    """Define previous VRS versions that are supported for computing digests and
    identifiers based on the current VRS model
    """

    V1_3 = "1.3"

    @classmethod
    def validate(cls, version):
        if version is not None and version not in cls.__members__.values():
            err_msg = f"Expected `PrevVrsVersion`, but got {version}"
            raise ValueError(err_msg)


ga4gh_compute_identifier_when = contextvars.ContextVar("ga4gh_compute_identifier_when")


class use_ga4gh_compute_identifier_when(ContextDecorator):
    """
    Context manager that defines when to compute identifiers
    for all operations within the context.  For example:

    with use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.GA4GH_INVALID):
        VCFAnnotator(...).annotate(...)

    Or:

    @use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.GA4GH_INVALID)
    def my_method():
    """

    def __init__(self, when: VrsObjectIdentifierIs):
        self.when = when
        self.token = None

    def __enter__(self):
        self.token = ga4gh_compute_identifier_when.set(self.when)

    def __exit__(self, exc_type, exc, exc_tb):
        ga4gh_compute_identifier_when.reset(self.token)


def is_ga4gh_identifier(ir):
    """

    >>> is_ga4gh_identifier("ga4gh:SQ.0123abcd")
    True

    >>> is_ga4gh_identifier("refseq:NM_01234.5")
    False

    >>> is_ga4gh_identifier(None)
    False

    """
    return str(get_pydantic_root(ir)).startswith(NS_W_SEP)


def ga4gh_identify(vro, in_place: str = 'default', as_version: PrevVrsVersion | None = None) -> str | None:
    """Return the GA4GH digest-based id for the object, as a CURIE
    (string).  Returns None if object is not identifiable.

    This function has three options for in_place editing of vro.id:
    - 'default': the standard identifier update behavior for GA4GH
        identifiable objects, this mode will update the vro.id
        field if the field is empty
    - 'always': this will update the vro.id field any time the
        identifier is computed (compute behavior is controlled by the
        use_ga4gh_compute_identifier_when context)
    - 'never': the vro.id field will not be edited in-place,
        even when empty

    If ``as_version`` is provided, other parameters are ignored and an identifier is
    returned following the conventions of the VRS version indicated by ``as_version_``.
    Raises ``ValueError`` if ``as_version`` is not a ``PrevVrsVersion``.

    >>> from ga4gh.core import ga4gh_identify
    >>> import ga4gh.vrs
    >>> location = ga4gh.vrs.models.SequenceLocation(start=44908821, end=44908822, sequenceReference=ga4gh.vrs.models.SequenceReference(refgetAccession="SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"))
    >>> ga4gh_identify(location)
    'ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT'
    """
    PrevVrsVersion.validate(as_version)

    if vro.is_ga4gh_identifiable():
        when_rule = ga4gh_compute_identifier_when.get(VrsObjectIdentifierIs.ANY)
        obj_id = None
        if when_rule == VrsObjectIdentifierIs.ANY:
            do_compute = True
        else:
            obj_id = getattr(vro, "id", None)
            if when_rule == VrsObjectIdentifierIs.MISSING:
                do_compute = obj_id is None or obj_id == ""
            else:  # VrsObjectIdentifierIs.GA4GH_INVALID
                do_compute = not vro.has_valid_ga4gh_id()

        if do_compute:
            obj_id = vro.get_or_create_ga4gh_identifier(in_place, as_version=as_version)

        return obj_id

    return None


def ga4gh_digest(vro: BaseModel, overwrite: bool = False, as_version: PrevVrsVersion | None = None) -> str | None:
    """Return the GA4GH digest for the object.

    Only GA4GH identifiable objects are GA4GH digestible.

    If ``as_version`` is provided, other parameters are ignored and a digest is returned
    following the conventions of the VRS version indicated by ``as_version_``.
    Raises ``ValueError`` if ``as_version`` is not a ``PrevVrsVersion``.

    >>> from ga4gh.core import ga4gh_digest
    >>> import ga4gh.vrs
    >>> location = ga4gh.vrs.models.SequenceLocation(start=44908821, end=44908822, sequenceReference=ga4gh.vrs.models.SequenceReference(refgetAccession="SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"))
    >>> ga4gh_digest(location)
    '4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT'
    """
    PrevVrsVersion.validate(as_version)

    if vro.is_ga4gh_identifiable():
        if as_version is None:
            return vro.get_or_create_digest(overwrite)
        else:
            return vro.compute_digest(as_version=as_version)
    else:
        return None


def ga4gh_serialize(obj: BaseModel, as_version: PrevVrsVersion | None = None) -> Optional[bytes]:
    """Serializes an object for use in computed digest computation.

    If ``as_version`` is provided, the returned serialization follows
    the conventions of the VRS version indicated by ``as_version_``.
    """
    PrevVrsVersion.validate(as_version)

    if as_version is None:
        return encode_canonical_json(obj.ga4gh_serialize())
    else:
        return obj.ga4gh_serialize_as_version(as_version)
