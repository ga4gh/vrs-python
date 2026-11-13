"""GA4GH VRS models

**This module should not be imported directly.**

Instead, users should use one of the following:

  * `from ga4gh.vrs import models`, and refer to models with the
    abbreviated name, e.g., `models.Allele` (recommended)

  * `import ga4gh.vrs`, and refer to models using the fully-qualified
    module name, e.g., `ga4gh.vrs.models.Allele`
"""

import inspect
import sys
from abc import ABC
from collections import OrderedDict
from enum import Enum
from types import UnionType
from typing import Annotated, Literal, get_args, get_origin

from canonicaljson import encode_canonical_json
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    RootModel,
    StringConstraints,
    ValidationInfo,
    field_validator,
)

from ga4gh.core import (
    CURIE_NAMESPACE,
    CURIE_SEP,
    GA4GH_IR_REGEXP,
    GA4GH_PREFIX_SEP,
    PrevVrsVersion,
    sha512t24u,
)
from ga4gh.core.models import (
    BaseModelForbidExtra,
    Element,
    Entity,
    iriReference,
)
from ga4gh.core.pydantic import get_pydantic_root, getattr_in


def flatten(vals):
    """Flatten vals recursively, lazily using yield"""

    def is_coll(thing) -> bool:
        """Return True if the thing looks like a collection.

        This is not exhaustive, do not use in general.
        """
        return type(thing) in [list, set]

    if is_coll(vals):
        for x in vals:
            yield from flatten(x)
    else:
        yield vals


def flatten_type(t):
    """Flatten a complex type into a list of constituent types."""
    origin = get_origin(t)
    args = get_args(t)

    if origin is Literal:
        return list(args)

    if origin in (UnionType, list):
        return list(flatten([flatten_type(sub_t) for sub_t in args]))
    return [t]


def overlaps(a: list, b: list) -> bool:
    """Return true if there are any elements in common between a and b"""
    return len(set(a).intersection(set(b))) > 0


def pydantic_class_refatt_map():
    """Build a map of class names to their field names that are referable types.

    As in, types with an identifier that can be referred to elsewhere, collapsed to that identifier and dereferenced.

    Returns a map like:

    {"Allele": ["location"], ...}
    """
    # Things defined here that are classes that inherit from BaseModel
    this_module = sys.modules[__name__]
    global_map = globals()
    model_classes = list(
        filter(
            lambda c: (
                inspect.isclass(c)
                and issubclass(c, BaseModel)
                and inspect.getmodule(c) == this_module
            ),
            [gl_name_value[1] for gl_name_value in global_map.items()],
        )
    )
    # Types directly reffable
    reffable_classes = list(
        filter(
            lambda c: (
                "id" in c.model_fields
                and getattr(c, "is_ga4gh_identifiable", lambda: False)()
            ),
            model_classes,
        )
    )
    # Types reffable because they are a union of reffable types
    union_reffable_classes = []
    for model_class in model_classes:
        if issubclass(model_class, RootModel):
            flattened_type_annotation = flatten_type(
                model_class.model_fields["root"].annotation
            )
            if overlaps(reffable_classes, flattened_type_annotation):
                union_reffable_classes.append(model_class)
    reffable_fields = {}
    # Find any field whose type is a subclass of a reffable type,
    # or which is a list that includes a reffable type.
    for model_class in model_classes:
        fields = model_class.model_fields
        class_reffable_fields = []
        for fieldname, field in fields.items():
            if fieldname == "root":
                continue
            field_type = field.annotation  # a typing or normal annotation like str
            # types can be raw class type annotation (int, str, dict, etc)
            # typing.Literal, typing.Union, typing.Optional
            # Use flatten_type to simplify these
            if any(
                [
                    rc in flatten_type(field_type)
                    for rc in (reffable_classes + union_reffable_classes)
                ]
            ):
                class_reffable_fields.append(fieldname)
        if len(class_reffable_fields) > 0:
            reffable_fields[model_class.__name__] = class_reffable_fields
    class_inherent = {}
    for model_class in model_classes:
        inherent = getattr_in(model_class, ["ga4gh", "inherent"])
        if inherent and len(inherent) > 0:
            class_inherent[model_class.__name__] = inherent
    return (reffable_classes, union_reffable_classes, reffable_fields, class_inherent)


class VrsType(str, Enum):
    """Define  VRS Types"""

    LEN_EXPR = "LengthExpression"
    REF_LEN_EXPR = "ReferenceLengthExpression"
    LIT_SEQ_EXPR = "LiteralSequenceExpression"
    SEQ_REF = "SequenceReference"
    SEQ_LOC = "SequenceLocation"
    ALLELE = "Allele"
    CIS_PHASED_BLOCK = "CisPhasedBlock"
    ADJACENCY = "Adjacency"
    TERMINUS = "Terminus"
    TRAVERSAL_BLOCK = "TraversalBlock"
    DERIVATIVE_MOL = "DerivativeMolecule"
    CN_COUNT = "CopyNumberCount"
    CN_CHANGE = "CopyNumberChange"


class Orientation(str, Enum):
    """The orientation of the molecular variation component."""

    FORWARD = "forward"
    REVERSE_COMPLEMENT = "reverse_complement"


class ResidueAlphabet(str, Enum):
    """The interpretation of the character codes referred to by the refget accession,
    where "aa" specifies an amino acid character set, and "na" specifies a nucleic acid
    character set.
    """

    AA = "aa"
    NA = "na"


class MoleculeType(str, Enum):
    """Molecule types as `defined by RefSeq <https://www.ncbi.nlm.nih.gov/books/NBK21091/>`_ (see Table 1)."""

    GENOMIC = "genomic"
    RNA = "RNA"
    MRNA = "mRNA"
    PROTEIN = "protein"


class CopyChange(str, Enum):
    """MUST use one of the defined enumerations that are based on the corresponding EFO
    ontological terms for copy number variation. See Implementation Guidance for more
    details.
    """

    COMPLETE_GENOMIC_LOSS = "complete genomic loss"
    HIGH_LEVEL_LOSS = "high-level loss"
    LOW_LEVEL_LOSS = "low-level loss"
    LOSS = "loss"
    REGIONAL_BASE_PLOIDY = "regional base ploidy"
    GAIN = "gain"
    LOW_LEVEL_GAIN = "low-level gain"
    HIGH_LEVEL_GAIN = "high-level gain"


class Syntax(str, Enum):
    """The syntax used to describe the variation. The value should be one of the
    supported syntaxes.
    """

    HGVS_C = "hgvs.c"
    HGVS_P = "hgvs.p"
    HGVS_G = "hgvs.g"
    HGVS_M = "hgvs.m"
    HGVS_N = "hgvs.n"
    HGVS_R = "hgvs.r"
    HGVS_ISCN = "iscn"
    GNOMAD = "gnomad"
    SPDI = "spdi"


def _recurse_ga4gh_serialize(obj, *, store: bool = True):
    """Serialize ``obj`` for GA4GH digest computation.

    Identifiable objects contribute a digest. Other value objects contribute
    their serialized inherent attributes. An identifiable object that already
    has ``digest`` set contributes that cached value.

    :param obj: Value to serialize
    :param store: When ``True`` (default), nested identifiable objects that do
        not yet have a digest store the computed digest on themselves. When
        ``False``, those digests are computed and returned without writing
        ``digest``.
    :return: A JSON-compatible value suitable for canonical encoding
    """
    if isinstance(obj, Ga4ghIdentifiableObject):
        return obj.get_or_create_digest(store=store)
    if isinstance(obj, _ValueObject):
        return obj.ga4gh_serialize(store=store)
    if isinstance(obj, RootModel):
        return _recurse_ga4gh_serialize(obj.model_dump(), store=store)
    if isinstance(obj, str):
        return obj
    if isinstance(obj, list):
        return [_recurse_ga4gh_serialize(x, store=store) for x in obj]
    return obj


class _ValueObject(Entity, ABC):
    """A contextual value whose equality is based on value, not identity.
    See https://en.wikipedia.org/wiki/Value_object for more on Value Objects.
    """

    def __hash__(self) -> int:
        return encode_canonical_json(self.ga4gh_serialize()).decode("utf-8").__hash__()

    def ga4gh_serialize(self, store: bool = True) -> dict:
        """Serialize the inherent attributes of this object for digest computation.

        :param store: When ``True`` (default), nested identifiable objects that
            do not yet have a digest store the computed digest on themselves.
            When ``False``, those digests are computed without writing
            ``digest``. An object that already has ``digest`` set contributes
            that cached value either way. ``recompute`` is not a serialization
            option: only ``get_or_create_digest`` decides whether this object's
            own cached digest is ignored, and nested cached digests are kept.
        :return: Mapping of inherent attribute names to serialized values
        """
        out = OrderedDict()
        for k in self.ga4gh.inherent:
            v = getattr(self, k)
            out[k] = _recurse_ga4gh_serialize(v, store=store)
        return out

    class ga4gh:  # noqa: N801
        inherent: list[str]

    @staticmethod
    def is_ga4gh_identifiable() -> bool:
        return False


class Ga4ghIdentifiableObject(_ValueObject, ABC):
    """A contextual value object for which a GA4GH computed identifier can be created.
    All GA4GH Identifiable Objects may have computed digests from the VRS Computed
    Identifier algorithm.
    """

    type: str
    digest: (
        Annotated[str, StringConstraints(pattern=r"^[0-9A-Za-z_\-]{32}$")] | None
    ) = Field(
        default=None,
        description="A sha512t24u digest created using the VRS Computed Identifier algorithm.",
    )

    def __lt__(self, other) -> bool:
        return self.get_or_create_digest() < other.get_or_create_digest()

    @staticmethod
    def is_ga4gh_identifiable() -> bool:
        return True

    def has_valid_ga4gh_id(self) -> bool | str | None:
        return self.id and GA4GH_IR_REGEXP.match(self.id) is not None

    def compute_digest(
        self, store: bool = True, as_version: PrevVrsVersion | None = None
    ) -> str:
        """Compute a sha512t24u digest using the VRS Computed Identifier algorithm.

        This always calculates a digest. Reuse of a cached ``digest`` is handled
        by ``get_or_create_digest``.

        :param store: When ``True`` (default), store the computed digest on this
            object, and let nested identifiable objects that lack a digest store
            theirs during serialization. When ``False``, return the digest
            without writing ``digest`` on this object or on nested identifiable
            objects. Nested objects that already have a digest still contribute
            that cached value.
        :param as_version: When provided, ``store`` is ignored and a digest is
            returned following the conventions of that VRS version.
        :raises AttributeError: If ``as_version`` is set and this class does not
            support prior version identifiers
        :return: The sha512t24u digest
        """
        if as_version is None:
            digest = sha512t24u(
                encode_canonical_json(self.ga4gh_serialize(store=store))
            )
            if store:
                self.digest = digest
        else:
            try:
                digest = sha512t24u(
                    self.ga4gh_serialize_as_version(as_version).encode("utf-8")
                )
            except AttributeError as e:
                msg = "This class does not support prior version identifiers."
                raise AttributeError(msg) from e
        return digest

    def get_or_create_ga4gh_identifier(
        self, in_place: str = "default", recompute: bool = False, as_version=None
    ) -> str:
        """Return a GA4GH computed identifier for this object.

        ``in_place="never"`` and ``in_place="always"`` return the same identifier.
        Both reuse this object's cached ``digest`` unless ``recompute`` is
        ``True``. ``never`` computes with ``store=False``, so ``id`` and every
        ``digest`` field stay as they were.

        :param in_place: Controls whether ``id`` is written on this object.
            ``default`` sets ``id`` only when it is empty. ``always`` sets
            ``id`` whenever the identifier is computed. ``never`` does not
            write ``id`` or any ``digest`` field.
        :param recompute: When ``True``, ignore this object's cached digest and
            compute a new one from current state. When ``False`` (default),
            return the identifier for ``digest`` unchanged when it is already
            set, including when other state has changed since that digest was
            stored. Nested identifiable objects that already have a digest keep
            that cached value; ``recompute`` does not refresh them.
        :param as_version: When provided, ``in_place`` and ``recompute`` are
            ignored and an identifier is returned following the conventions of
            that VRS version.
        :raises ValueError: If ``in_place`` is not ``default``, ``always``, or
            ``never``
        :raises AttributeError: If ``as_version`` is set and this class does not
            support prior version identifiers
        :return: The GA4GH computed identifier CURIE
        """
        if as_version is not None:
            return self.compute_ga4gh_identifier(as_version=as_version)

        if in_place == "default":
            if self.id is None:
                self.id = self.compute_ga4gh_identifier(recompute)
        elif in_place == "always":
            self.id = self.compute_ga4gh_identifier(recompute)
        elif in_place == "never":
            return self.compute_ga4gh_identifier(recompute, store=False)
        else:
            msg = "Expected 'in_place' to be one of 'default', 'always', or 'never'"
            raise ValueError(msg)

        if self.has_valid_ga4gh_id():
            return self.id
        else:
            return self.compute_ga4gh_identifier(recompute)

    def compute_ga4gh_identifier(
        self, recompute: bool = False, as_version=None, store: bool = True
    ) -> str:
        """Return a GA4GH computed identifier.

        :param recompute: When ``True``, ignore this object's cached digest and
            compute a new one. When ``False`` (default), reuse ``digest`` if it
            is already set. Nested identifiable objects that already have a
            digest keep that cached value.
        :param as_version: When provided, ``recompute`` and ``store`` are
            ignored and an identifier is returned following the conventions of
            that VRS version.
        :param store: When ``True`` (default), persist a newly computed digest
            on this object and on nested identifiable objects that do not yet
            have one. When ``False``, compute the identifier without writing
            ``digest`` fields. The returned identifier matches the ``store=True``
            result for the same ``recompute`` value.
        :raises AttributeError: If ``as_version`` is set and this class does not
            support prior version identifiers
        :return: The GA4GH computed identifier CURIE
        """
        if as_version is None:
            digest = self.get_or_create_digest(recompute, store=store)
            return f"{CURIE_NAMESPACE}{CURIE_SEP}{self.ga4gh.prefix}{GA4GH_PREFIX_SEP}{digest}"
        else:
            digest = self.compute_digest(as_version=as_version)
            return f"{CURIE_NAMESPACE}{CURIE_SEP}{self.ga4gh.priorPrefix[as_version]}{GA4GH_PREFIX_SEP}{digest}"

    def get_or_create_digest(self, recompute: bool = False, store: bool = True) -> str:
        """Return this object's sha512t24u digest, computing it when needed.

        :param recompute: When ``True``, ignore this object's cached digest and
            compute a new one from current state. When ``False`` (default),
            return ``digest`` as stored if it is set. This does not recompute
            nested digests that are already cached. A cached digest returned
            on this path is not written again.
        :param store: When ``True`` (default), a newly computed digest is
            written to this object, and nested identifiable objects that lack a
            digest store theirs. When ``False``, the digest is computed and
            returned without writing ``digest`` on this object or on nested
            objects.
        :return: The sha512t24u digest
        """
        if self.digest is None or recompute:
            return self.compute_digest(store=store)
        return self.digest

    class ga4gh(_ValueObject.ga4gh):  # noqa: N801
        prefix: str


class Expression(Element, BaseModelForbidExtra):
    """Representation of a variation by a specified nomenclature or syntax for a
    Variation object. Common examples of expressions for the description of molecular
    variation include the HGVS and ISCN nomenclatures.
    """

    model_config = ConfigDict(use_enum_values=True)

    syntax: Syntax = Field(
        ...,
        description="The syntax used to describe the variation. The value should be one of the supported syntaxes.",
    )
    value: str = Field(
        ...,
        description="The expression of the variation in the specified syntax. The value should be a valid expression in the specified syntax.",
    )
    syntax_version: str | None = Field(
        default=None,
        description="The version of the syntax used to describe the variation. This is particularly important for HGVS expressions, as the syntax has evolved over time.",
    )


#########################################
# vrs numerics, comparators, and ranges
#########################################


class Range(RootModel):
    """An inclusive range of values bounded by one or more integers."""

    root: list[int | None] = Field(
        ...,
        json_schema_extra={
            "description": "An inclusive range of values bounded by one or more integers."
        },
        max_length=2,
        min_length=2,
    )

    @field_validator("root", mode="after")
    def validate_range(cls, v: list[int | None]) -> list[int | None]:  # noqa: N805
        """Validate range values

        :param v: Root value
        :raises ValueError: If ``root`` does not include at least one integer or if
            the first element in ``root`` is greater than the second element in ``root``
        :return: Inclusive range
        """
        if v.count(None) == 2:
            err_msg = "Must provide at least one integer."
            raise ValueError(err_msg)

        if (v[0] is not None and v[1] is not None) and (v[0] > v[1]):
            err_msg = (
                "The first integer must be less than or equal to the second integer."
            )
            raise ValueError(err_msg)

        return v


class residue(RootModel):
    """A character representing a specific residue (i.e., molecular species) or
    groupings of these ("ambiguity codes"), using `one-letter IUPAC abbreviations
    <https://en.wikipedia.org/wiki/International_Union_of_Pure_and_Applied_Chemistry#Amino_acid_and_nucleotide_base_codes>`_
    for nucleic acids and amino acids.
    """

    root: Annotated[str, StringConstraints(pattern=r"[A-Z*\-]")] = Field(
        ...,
        json_schema_extra={
            "description": 'A character representing a specific residue (i.e., molecular species) or groupings of these ("ambiguity codes"), using [one-letter IUPAC abbreviations](https://en.wikipedia.org/wiki/International_Union_of_Pure_and_Applied_Chemistry#Amino_acid_and_nucleotide_base_codes) for nucleic acids and amino acids.'
        },
    )


class sequenceString(RootModel):
    """A character string of `Residues` that represents a biological sequence using the
    conventional sequence order (5'-to-3' for nucleic acid sequences, and
    amino-to-carboxyl for amino acid sequences). IUPAC ambiguity codes are permitted in
    Sequence Strings.
    """

    root: Annotated[str, StringConstraints(pattern=r"^[A-Z*\-]*$")] = Field(
        ...,
        json_schema_extra={
            "description": "A character string of Residues that represents a biological sequence using the conventional sequence order (5'-to-3' for nucleic acid sequences, and amino-to-carboxyl for amino acid sequences). IUPAC ambiguity codes are permitted in Sequence Strings."
        },
    )


#########################################
# vrs sequence expression
#########################################


class LengthExpression(_ValueObject, BaseModelForbidExtra):
    """A sequence expressed only by its length."""

    type: Literal["LengthExpression"] = Field(
        default=VrsType.LEN_EXPR.value,
        description=f'MUST be "{VrsType.LEN_EXPR.value}"',
    )
    length: Range | int | None = Field(
        default=None, description="The length of the sequence."
    )

    class ga4gh(_ValueObject.ga4gh):
        inherent = ["length", "type"]


class ReferenceLengthExpression(_ValueObject, BaseModelForbidExtra):
    """An expression of a length of a sequence from a repeating reference."""

    type: Literal["ReferenceLengthExpression"] = Field(
        default=VrsType.REF_LEN_EXPR.value,
        description=f'MUST be "{VrsType.REF_LEN_EXPR.value}"',
    )
    length: Range | int = Field(
        ..., description="The number of residues in the expressed sequence."
    )
    sequence: sequenceString | None = Field(
        default=None,
        description="the literal Sequence encoded by the Reference Length Expression.",
    )
    repeatSubunitLength: int = Field(
        ..., description="The number of residues in the repeat subunit."
    )

    class ga4gh(_ValueObject.ga4gh):
        inherent = ["length", "repeatSubunitLength", "type"]


class LiteralSequenceExpression(_ValueObject, BaseModelForbidExtra):
    """An explicit expression of a Sequence."""

    type: Literal["LiteralSequenceExpression"] = Field(
        default=VrsType.LIT_SEQ_EXPR.value,
        description=f'MUST be "{VrsType.LIT_SEQ_EXPR.value}"',
    )
    sequence: sequenceString = Field(..., description="the literal sequence")

    class ga4gh(_ValueObject.ga4gh):
        inherent = ["sequence", "type"]


#########################################
# vrs location
#########################################


class SequenceReference(_ValueObject, BaseModelForbidExtra):
    """A sequence of nucleic or amino acid character codes."""

    model_config = ConfigDict(use_enum_values=True)

    type: Literal["SequenceReference"] = Field(
        default=VrsType.SEQ_REF.value, description=f'MUST be "{VrsType.SEQ_REF.value}"'
    )
    refgetAccession: Annotated[
        str, StringConstraints(pattern=r"^SQ.[0-9A-Za-z_\-]{32}$")
    ] = Field(
        ...,
        description="A [GA4GH RefGet](http://samtools.github.io/hts-specs/refget.html) identifier for the referenced sequence, using the sha512t24u digest.",
    )
    residueAlphabet: ResidueAlphabet | None = Field(
        default=None,
        description='The interpretation of the character codes referred to by the refget accession, where "aa" specifies an amino acid character set, and "na" specifies a nucleic acid character set.',
    )
    circular: bool | None = Field(
        default=None,
        description="A boolean indicating whether the molecule represented by the sequence is circular (true) or linear (false).",
    )
    sequence: sequenceString | None = Field(
        default=None,
        description="A sequenceString that is a literal representation of the referenced sequence.",
    )
    moleculeType: MoleculeType | None = Field(
        default=None,
        description="Molecule types as [defined by RefSeq](https://www.ncbi.nlm.nih.gov/books/NBK21091/) (see Table 1). MUST be one of 'genomic', 'RNA', 'mRNA', or 'protein'.",
    )

    class ga4gh(_ValueObject.ga4gh):
        inherent = ["refgetAccession", "type"]


class SequenceLocation(Ga4ghIdentifiableObject, BaseModelForbidExtra):
    """A `Location` defined by an interval on a `Sequence`."""

    type: Literal["SequenceLocation"] = Field(
        default=VrsType.SEQ_LOC.value, description=f'MUST be "{VrsType.SEQ_LOC.value}"'
    )
    sequenceReference: iriReference | SequenceReference | None = Field(
        default=None,
        description="A reference to a SequenceReference on which the location is defined.",
    )
    start: Range | int | None = Field(
        default=None,
        description="The start coordinate or range of the SequenceLocation. The minimum value of this coordinate or range is 0. For locations on linear sequences, this MUST represent a coordinate or range less than or equal to the value of `end`. For circular sequences, `start` is greater than `end` when the location spans the sequence 0 coordinate.",
    )
    end: Range | int | None = Field(
        default=None,
        description="The end coordinate or range of the SequenceLocation. The minimum value of this coordinate or range is 0. For locations on linear sequences, this MUST represent a coordinate or range greater than or equal to the value of `start`. For circular sequences, `end` is less than `start` when the location spans the sequence 0 coordinate.",
    )
    sequence: sequenceString | None = Field(
        default=None,
        description="The literal sequence encoded by the `sequenceReference` at these coordinates.",
    )

    @field_validator("start", "end", mode="after")
    def validate_start_end(
        cls, v: Range | int | None, info: ValidationInfo
    ) -> Range | int | None:
        """Validate ``start`` and ``end`` fields

        :param v: ``start`` or ``end`` value
        :param info: Validation info
        :raises ValueError: If ``start`` or ``end`` has a value less than 0
        :return: Sequence Location
        """
        if v is not None:
            if isinstance(v, int):
                int_values = [v]
            else:
                int_values = [val for val in v.root if val is not None]

            if any(int_val < 0 for int_val in int_values):
                err_msg = f"The minimum value of `{info.field_name}` is 0."
                raise ValueError(err_msg)
        return v

    def ga4gh_serialize_as_version(self, as_version: PrevVrsVersion) -> str:
        """Return a serialized string following the conventions for SequenceLocation
        serialization as defined in the VRS version specified by ``as_version``.

        :raises ValueError: If ``sequenceReference`` is not a ``SequenceReference``
            object; ``start`` or ``end`` are not an int or list.
        """
        if as_version == PrevVrsVersion.V1_3:
            if not isinstance(self.sequenceReference, SequenceReference):
                err_msg = "Must provide `sequenceReference` and it must be a valid `SequenceReference`"
                raise ValueError(err_msg)

            out = []
            for value in [self.start, self.end]:
                value = get_pydantic_root(value)
                if isinstance(value, int):
                    result = f'{{"type":"Number","value":{value}}}'
                elif isinstance(value, list):
                    if value[0] is None:
                        result = f'{{"comparator":"<=","type":"IndefiniteRange","value":{value[1]}}}'
                    elif value[1] is None:
                        result = f'{{"comparator":">=","type":"IndefiniteRange","value":{value[0]}}}'
                    else:
                        result = f'{{"max":{value[1]},"min":{value[0]},"type":"DefiniteRange"}}'
                else:
                    msg = f"{value} is not int or list."
                    raise ValueError(msg)
                out.append(result)
            return f'{{"interval":{{"end":{out[1]},"start":{out[0]},"type":"SequenceInterval"}},"sequence_id":"{self.sequenceReference.refgetAccession.split(".")[1]}","type":"SequenceLocation"}}'
        msg = f"Received an unexpected value for `as_version`: {as_version}. MUST be an instance of `PrevVrsVersion`."
        raise TypeError(msg)

    def get_refget_accession(self) -> str | None:
        if isinstance(self.sequenceReference, SequenceReference):
            return self.sequenceReference.refgetAccession
        if isinstance(self.sequenceReference, iriReference):
            return self.sequenceReference.root
        return None

    class ga4gh(Ga4ghIdentifiableObject.ga4gh):  # noqa: N801
        prefix = "SL"
        priorPrefix = {PrevVrsVersion.V1_3.value: "VSL"}  # noqa: N815
        inherent = ["end", "sequenceReference", "start", "type"]


#########################################
# base variation
#########################################


class _VariationBase(Ga4ghIdentifiableObject, ABC):
    """Base class for variation"""

    expressions: list[Expression] | None = None


#########################################
# vrs molecular variation
#########################################


class Allele(_VariationBase, BaseModelForbidExtra):
    """The state of a molecule at a `Location`."""

    type: Literal["Allele"] = Field(
        default=VrsType.ALLELE.value, description=f'MUST be "{VrsType.ALLELE.value}"'
    )
    location: iriReference | SequenceLocation = Field(
        ..., description="The location of the Allele"
    )
    state: LiteralSequenceExpression | ReferenceLengthExpression | LengthExpression = (
        Field(..., description="An expression of the sequence state")
    )

    def ga4gh_serialize_as_version(self, as_version: PrevVrsVersion) -> str:
        """Return a serialized string following the conventions for
        Allele serialization as defined in the VRS version specified by 'as_version`.

        :raises ValueError: If ``state`` is not a ``LiteralSequenceExpression`` or
            ``ReferenceLengthExpression``; ``state.sequence`` is null.
        """
        location_digest = self.location.compute_digest(as_version=as_version)

        if not isinstance(
            self.state, LiteralSequenceExpression | ReferenceLengthExpression
        ):
            err_msg = "Only `LiteralSequenceExpression` and `ReferenceLengthExpression` are supported for previous versions of VRS"
            raise ValueError(err_msg)

        sequence = get_pydantic_root(self.state.sequence)

        if sequence is None:
            raise ValueError("State sequence attribute must be defined.")

        if as_version == PrevVrsVersion.V1_3:
            return f'{{"location":"{location_digest}","state":{{"sequence":"{sequence}","type":"LiteralSequenceExpression"}},"type":"Allele"}}'
        msg = f"Received an unexpected value for `as_version`: {as_version}. MUST be an instance of `PrevVrsVersion`."
        raise TypeError(msg)

    class ga4gh(Ga4ghIdentifiableObject.ga4gh):  # noqa: N801
        prefix = "VA"
        priorPrefix = {PrevVrsVersion.V1_3.value: "VA"}  # noqa: N815
        inherent = ["location", "state", "type"]


class CisPhasedBlock(_VariationBase, BaseModelForbidExtra):
    """An ordered set of co-occurring `Variation` on the same molecule."""

    type: Literal["CisPhasedBlock"] = Field(
        default=VrsType.CIS_PHASED_BLOCK.value,
        description=f'MUST be "{VrsType.CIS_PHASED_BLOCK.value}"',
    )
    members: list[Allele | iriReference] = Field(
        ...,
        description="A list of Alleles that are found in-cis on a shared molecule.",
        min_length=2,
    )
    sequenceReference: SequenceReference | None = Field(
        default=None,
        description="An optional Sequence Reference on which all of the in-cis Alleles are found. When defined, this may be used to implicitly define the `sequenceReference` attribute for each of the CisPhasedBlock member Alleles.",
    )

    def ga4gh_serialize(self, store: bool = True) -> dict:
        """Serialize a CisPhasedBlock, sorting member digests.

        Member order is not inherent. Digests are sorted so two blocks with
        the same members produce the same identifier. ``store`` is forwarded
        to the shared serializer so ``store=False`` still sorts, and still
        avoids writing nested digests.

        :param store: When ``True`` (default), nested identifiable objects that
            do not yet have a digest store the computed digest on themselves.
            When ``False``, those digests are computed without writing
            ``digest``. Cached digests are reused either way.
        :return: Serialized inherent attributes, with ``members`` sorted by digest
        """
        out = _ValueObject.ga4gh_serialize(self, store=store)
        out["members"] = sorted(out["members"])
        return out

    class ga4gh(Ga4ghIdentifiableObject.ga4gh):
        prefix = "CPB"
        inherent = ["members", "type"]


#########################################
# vrs structural variation (under active discussion)
#########################################


class Adjacency(_VariationBase, BaseModelForbidExtra):
    """The `Adjacency` class represents the adjoining of the end of a sequence with the
    beginning of an adjacent sequence, potentially with an intervening linker sequence.
    """

    type: Literal["Adjacency"] = Field(
        default=VrsType.ADJACENCY.value,
        description=f'MUST be "{VrsType.ADJACENCY.value}".',
    )
    adjoinedSequences: list[iriReference | SequenceLocation] = Field(
        ...,
        description="The terminal sequence or pair of adjoined sequences that defines in the adjacency.",
        min_length=2,
        max_length=2,
    )
    linker: (
        LiteralSequenceExpression | ReferenceLengthExpression | LengthExpression | None
    ) = Field(
        default=None, description="The sequence found between adjoined sequences."
    )
    homology: bool | None = Field(
        default=None,
        description="A flag indicating if coordinate ambiguity in the adjoined sequences is from sequence homology (true) or other uncertainty, such as instrument ambiguity (false).",
    )

    @field_validator("adjoinedSequences", mode="after")
    def validate_adjoined_sequences(cls, v) -> list[iriReference | SequenceLocation]:
        """Ensure ``adjoinedSequences`` do not have both ``start`` and ``end``

        :raises ValueError: If an adjoined sequence has both ``start`` and ``end``
        :return: Adjoined sequences represented as iri reference or sequence location
        """
        for adjoined_seq in v:
            if isinstance(adjoined_seq, SequenceLocation):
                if adjoined_seq.start and adjoined_seq.end:
                    err_msg = "Adjoined sequence must not have both `start` and `end`."
                    raise ValueError(err_msg)
        return v

    class ga4gh(Ga4ghIdentifiableObject.ga4gh):
        prefix = "AJ"
        inherent = ["adjoinedSequences", "linker", "type"]


class Terminus(_VariationBase, BaseModelForbidExtra):
    """The `Terminus` data class provides a structure for describing the end
    (terminus) of a sequence. Structurally similar to Adjacency but the linker sequence
    is not allowed and it removes the unnecessary array structure.
    """

    type: Literal["Terminus"] = Field(
        default=VrsType.TERMINUS.value,
        description=f'MUST be "{VrsType.TERMINUS.value}".',
    )
    location: iriReference | SequenceLocation = Field(
        ..., description="The location of the terminus."
    )

    class ga4gh(Ga4ghIdentifiableObject.ga4gh):  # noqa: N815
        prefix = "TM"
        inherent = ["location", "type"]


class TraversalBlock(_ValueObject, BaseModelForbidExtra):
    """A component used to describe the orientation of applicable molecular variation
    within a DerivativeMolecule.
    """

    model_config = ConfigDict(use_enum_values=True)

    type: Literal["TraversalBlock"] = Field(
        default=VrsType.TRAVERSAL_BLOCK.value,
        description=f'MUST be "{VrsType.TRAVERSAL_BLOCK.value}".',
    )
    orientation: Orientation | None = Field(
        default=None,
        description="The orientation of the molecular variation component.",
    )

    component: Adjacency | None = Field(
        default=None, description="The unoriented molecular variation component."
    )

    class ga4gh(_ValueObject.ga4gh):
        inherent = ["component", "orientation", "type"]


class DerivativeMolecule(_VariationBase, BaseModelForbidExtra):
    """The "Derivative Molecule" data class is a structure for describing a derivate
    molecule composed from multiple sequence components.
    """

    type: Literal["DerivativeMolecule"] = Field(
        default=VrsType.DERIVATIVE_MOL.value,
        description=f'MUST be "{VrsType.DERIVATIVE_MOL.value}".',
    )
    components: list[
        iriReference | Allele | CisPhasedBlock | Terminus | TraversalBlock
    ] = Field(
        ...,
        description="The molecular components that constitute the derivative molecule.",
        min_length=2,
    )
    circular: bool | None = Field(
        default=None,
        description="A boolean indicating whether the molecule represented by the sequence is circular (true) or linear (false).",
    )

    class ga4gh(Ga4ghIdentifiableObject.ga4gh):  # noqa: N815
        prefix = "DM"
        inherent = ["components", "type"]


#########################################
# vrs systemic variation
#########################################


class CopyNumberCount(_VariationBase, BaseModelForbidExtra):
    """The absolute count of discrete copies of a `Location`, within a system
    (e.g. genome, cell, etc.).
    """

    type: Literal["CopyNumberCount"] = Field(
        default=VrsType.CN_COUNT.value,
        description=f'MUST be "{VrsType.CN_COUNT.value}"',
    )
    location: iriReference | SequenceLocation = Field(
        ...,
        description="The location of the subject of the copy count.",
    )
    copies: Range | int = Field(
        ..., description="The integral number of copies of the subject in a system"
    )

    class ga4gh(Ga4ghIdentifiableObject.ga4gh):  # noqa: N815
        prefix = "CN"
        inherent = ["copies", "location", "type"]


class CopyNumberChange(_VariationBase, BaseModelForbidExtra):
    """An assessment of the copy number of a `Location` within a system
    (e.g. genome, cell, etc.) relative to a baseline ploidy.
    """

    model_config = ConfigDict(use_enum_values=True)

    type: Literal["CopyNumberChange"] = Field(
        default=VrsType.CN_CHANGE.value,
        description=f'MUST be "{VrsType.CN_CHANGE.value}"',
    )
    location: iriReference | SequenceLocation = Field(
        ...,
        description="The location of the subject of the copy change.",
    )
    copyChange: CopyChange = Field(
        ...,
        description="MUST use one of the defined enumerations that are based on the corresponding EFO ontological terms for copy number variation. See Implementation Guidance for more details.",
    )

    class ga4gh(Ga4ghIdentifiableObject.ga4gh):
        prefix = "CX"
        inherent = ["copyChange", "location", "type"]


#########################################
# vrs kinds of variation, expression, and location
#########################################


class MolecularVariation(RootModel):
    """A `variation` on a contiguous molecule."""

    root: Allele | CisPhasedBlock | Adjacency | Terminus | DerivativeMolecule = Field(
        ...,
        json_schema_extra={"description": "A `variation` on a contiguous molecule."},
        discriminator="type",
    )


class SequenceExpression(RootModel):
    """An expression describing a `Sequence`."""

    root: LiteralSequenceExpression | ReferenceLengthExpression | LengthExpression = (
        Field(
            ...,
            json_schema_extra={"description": "An expression describing a `Sequence`."},
            discriminator="type",
        )
    )


class Location(RootModel):
    """A contiguous segment of a biological sequence."""

    root: SequenceLocation = Field(
        ...,
        json_schema_extra={
            "description": "A contiguous segment of a biological sequence."
        },
        discriminator="type",
    )


class Variation(RootModel):
    """A representation of the state of one or more biomolecules."""

    root: (
        Allele
        | CisPhasedBlock
        | Adjacency
        | Terminus
        | DerivativeMolecule
        | CopyNumberChange
        | CopyNumberCount
    ) = Field(
        ...,
        json_schema_extra={
            "description": "A representation of the state of one or more biomolecules."
        },
        discriminator="type",
    )


class SystemicVariation(RootModel):
    """A Variation of multiple molecules in the context of a system, e.g. a genome,
    sample, or homologous chromosomes.
    """

    root: CopyNumberChange | CopyNumberCount = Field(
        ...,
        json_schema_extra={
            "description": "A Variation of multiple molecules in the context of a system, e.g. a genome, sample, or homologous chromosomes."
        },
        discriminator="type",
    )


# At end so classes exist
(reffable_classes, union_reffable_classes, class_refatt_map, class_inherent) = (
    pydantic_class_refatt_map()
)
