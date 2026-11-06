"""Serialize, digest, and identify GA4GH models.

For an Allele without cached IDs or digests, the main steps are
(intermediate helper calls omitted)::

    ga4gh_identify(allele)
    + allele.compute_ga4gh_identifier()
      + allele.compute_digest()
        + allele.ga4gh_serialize()
          + location.get_or_create_digest()  # identifiable child
          + state.ga4gh_serialize()          # inline value object
        + encode_canonical_json(...)
        + sha512t24u(...)

Given an Allele, call the public helpers::

    from ga4gh.core import ga4gh_digest, ga4gh_identify, ga4gh_serialize

    identifier = ga4gh_identify(allele, in_place="never")  # GA4GH CURIE
    digest = ga4gh_digest(allele)  # digest component
    serialized = ga4gh_serialize(allele)  # canonical JSON bytes

``in_place="never"`` leaves the ID unchanged; ``"always"`` stores the result;
``"default"`` stores it only when the ID is None. Digests may still be cached.

To reuse any nonempty ID::

    from ga4gh.core import VrsObjectIdentifierIs, use_ga4gh_compute_identifier_when

    with use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.MISSING):
        identifier = ga4gh_identify(allele)

``GA4GH_INVALID`` instead requires an ID matching the GA4GH pattern. ``ANY``
always calls the model identifier method, which may reuse cached IDs or digests.
After editing an allele, refresh changed identifiable children first, then::

    identifier = allele.get_or_create_ga4gh_identifier(
        in_place="always", recompute=True
    )
"""

import contextvars
import re
from contextlib import ContextDecorator
from enum import Enum, IntEnum
from types import TracebackType

from canonicaljson import encode_canonical_json
from pydantic import BaseModel

from ga4gh.core.pydantic import get_pydantic_root

__all__ = ["ga4gh_digest", "ga4gh_identify", "ga4gh_serialize", "is_ga4gh_identifier"]

CURIE_NAMESPACE = "ga4gh"
CURIE_SEP = ":"
GA4GH_PREFIX_SEP = "."

GA4GH_IR_REGEXP = re.compile(r"^ga4gh:(?P<type>[^.]+)\.(?P<digest>[0-9A-Za-z_\-]{32})$")
GA4GH_DIGEST_REGEXP = re.compile(r"^[0-9A-Za-z_\-]{32}$")

NS_W_SEP = f"{CURIE_NAMESPACE}{CURIE_SEP}"


class VrsObjectIdentifierIs(IntEnum):
    """Control when ga4gh_identify calls the model identifier method.

    * ``ANY``: Always call the method (default).
    * ``GA4GH_INVALID``: Call it when the ID is missing or syntactically invalid.
    * ``MISSING``: Call it only when the ID is None or empty.

    The default policy calls the model identifier method. That method may reuse an
    existing valid identifier or digest; this policy does not force recomputation. Where
    the source of inputs to `ga4gh_identify` are well controlled, for example when
    annotating a VCF file with VRS IDs, using `MISSING` can improve performance.
    """

    ANY = 0
    GA4GH_INVALID = 1
    MISSING = 2


class PrevVrsVersion(str, Enum):
    """Define prior VRS versions supported for digest and identifier computation."""

    V1_3 = "1.3"

    @classmethod
    def validate(cls, version: object) -> None:
        """Validate a supported prior VRS version.

        :param version: Version value to validate; None selects the current version.
        :raises ValueError: If the version is unsupported.
        """
        if version is not None and version not in cls.__members__.values():
            err_msg = f"Expected `PrevVrsVersion`, but got {version}"
            raise ValueError(err_msg)


ga4gh_compute_identifier_when = contextvars.ContextVar("ga4gh_compute_identifier_when")


class use_ga4gh_compute_identifier_when(ContextDecorator):  # noqa: N801
    """Apply an identifier policy within a context or decorated function.

    As a context manager::

        with use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.GA4GH_INVALID):
            VCFAnnotator(...).annotate(...)

    As a decorator::

        @use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.GA4GH_INVALID)
        def my_method(): ...
    """

    def __init__(self, when: VrsObjectIdentifierIs) -> None:
        """Set the identifier computation policy.

        :param when: Identifier computation policy for the context.
        """
        self.when = when
        self.token = None

    def __enter__(self) -> None:
        """Apply the identifier computation policy."""
        self.token = ga4gh_compute_identifier_when.set(self.when)

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        """Restore the previous identifier computation policy.

        :param exc_type: Exception class from the context body, if any.
        :param exc: Exception instance from the context body, if any.
        :param exc_tb: Traceback from the context body, if any.
        """
        ga4gh_compute_identifier_when.reset(self.token)


def is_ga4gh_identifier(ir: object) -> bool:
    """Check whether a value starts with the GA4GH namespace prefix

    >>> is_ga4gh_identifier("ga4gh:SQ.0123abcd")
    True

    >>> is_ga4gh_identifier("refseq:NM_01234.5")
    False

    >>> is_ga4gh_identifier(None)
    False

    :param ir: Value or root model to inspect.
    :return: Whether the value has the GA4GH namespace prefix.
    """
    return str(get_pydantic_root(ir)).startswith(NS_W_SEP)


def ga4gh_identify(
    vro: BaseModel,
    in_place: str = "default",
    as_version: PrevVrsVersion | None = None,
) -> str | None:
    """Return a GA4GH identifier, or None for a non-identifiable object.

    This function has three options for in_place editing of vro.id:

    - 'default': the standard identifier update behavior for GA4GH
        identifiable objects, this mode will update the vro.id
        field only when it is None
    - 'always': this will update the vro.id field any time the
        identifier is computed (compute behavior is controlled by the
        use_ga4gh_compute_identifier_when context)
    - 'never': the vro.id field will not be edited in-place,
        even when empty

    The context policy may return an existing ID without computing one. When computation
    occurs, ``as_version`` selects legacy serialization and bypasses in-place ID
    updates; it does not override the context policy. Accepts a supported PrevVrsVersion
    member or its string value. Raises ValueError for other version values.

    >>> from ga4gh.core import ga4gh_identify
    >>> import ga4gh.vrs
    >>> location = ga4gh.vrs.models.SequenceLocation(
    ...     start=44908821,
    ...     end=44908822,
    ...     sequenceReference=ga4gh.vrs.models.SequenceReference(
    ...         refgetAccession="SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
    ...     ),
    ... )
    >>> ga4gh_identify(location)
    'ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT'

    :param vro: GA4GH model implementing is_ga4gh_identifiable and the digest
        and identifier methods used for identifiable objects.
    :param in_place: Identifier update policy: default, always, or never.
    :param as_version: Optional prior VRS version to use.
    :return: Identifier, or None for a non-identifiable model.
    :raises ValueError: If as_version is unsupported.
    :raises AttributeError: If the model lacks the required GA4GH methods.
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


def ga4gh_digest(
    vro: BaseModel, overwrite: bool = False, as_version: PrevVrsVersion | None = None
) -> str | None:
    """Return the GA4GH digest for the object.

    Only GA4GH identifiable objects are GA4GH digestible.

    If ``as_version`` is provided, other parameters are ignored and a digest is returned
    following the conventions of the VRS version indicated by ``as_version``. Accepts a
    supported PrevVrsVersion member or its string value. Raises ValueError for other
    version values.

    >>> from ga4gh.core import ga4gh_digest
    >>> import ga4gh.vrs
    >>> location = ga4gh.vrs.models.SequenceLocation(
    ...     start=44908821,
    ...     end=44908822,
    ...     sequenceReference=ga4gh.vrs.models.SequenceReference(
    ...         refgetAccession="SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
    ...     ),
    ... )
    >>> ga4gh_digest(location)
    '4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT'

    :param vro: GA4GH model implementing is_ga4gh_identifiable and the digest
        and identifier methods used for identifiable objects.
    :param overwrite: Recompute an existing digest when true.
    :param as_version: Optional prior VRS version to use.
    :return: Digest, or None for a non-identifiable model.
    :raises ValueError: If as_version is unsupported.
    :raises AttributeError: If the model lacks the required GA4GH methods.
    """
    PrevVrsVersion.validate(as_version)

    if vro.is_ga4gh_identifiable():
        if as_version is None:
            return vro.get_or_create_digest(overwrite)
        return vro.compute_digest(as_version=as_version)
    return None


def ga4gh_serialize(
    obj: BaseModel, as_version: PrevVrsVersion | None = None
) -> bytes | str:
    """Serialize an object for use in computed digest computation.

    If ``as_version`` is provided, the returned serialization follows the conventions of
    the VRS version indicated by ``as_version``.

    :param obj: Model implementing ga4gh_serialize, and ga4gh_serialize_as_version
        when legacy serialization is requested.
    :param as_version: Optional prior VRS version to use.
    :return: Canonical JSON bytes, or legacy serialization text.
    :raises ValueError: If as_version is unsupported.
    :raises AttributeError: If the model lacks the required GA4GH methods.
    """
    PrevVrsVersion.validate(as_version)

    if as_version is None:
        return encode_canonical_json(obj.ga4gh_serialize())
    return obj.ga4gh_serialize_as_version(as_version)
