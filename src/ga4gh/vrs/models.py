"""GA4GH VRS models

**This module should not be imported directly.**

Instead, users should use one of the following:

  * `from ga4gh.vrs import models`, and refer to models with the
    abbreviated name, e.g., `models.Allele` (recommended)

  * `import ga4gh.vrs`, and refer to models using the fully-qualified
    module name, e.g., `ga4gh.vrs.models.Allele`
"""
from abc import ABC
from typing import List, Literal, Optional, Union, Dict, Annotated
from collections import OrderedDict
from enum import Enum
import inspect
import sys
from ga4gh.core import (
    sha512t24u,
    GA4GH_PREFIX_SEP,
    CURIE_SEP,
    CURIE_NAMESPACE,
    GA4GH_IR_REGEXP,
    PrevVrsVersion
)
from ga4gh.core.pydantic import get_pydantic_root

from canonicaljson import encode_canonical_json
from pydantic import BaseModel, Field, RootModel, StringConstraints, ConfigDict

from ga4gh.core.pydantic import (
    getattr_in
)
from ga4gh.core.entity_models import IRI, Expression, DomainEntity


def flatten(vals):
    """
    Flattens vals recursively, lazily using yield
    """
    def is_coll(thing):
        """
        Return True if the thing looks like a collection.
        This is not exhaustive, do not use in general.
        """
        return type(thing) in [list, set]
    if is_coll(vals):
        for x in vals:
            for fx in flatten(x):
                yield fx
    else:
        yield vals


def flatten_type(t):
    """
    Flattens a complex type into a list of constituent types.
    """
    if hasattr(t, '__dict__') and '__origin__' in t.__dict__:
        if t.__origin__ == Literal:
            return list(t.__args__)
        elif (t.__origin__ == Union
              or issubclass(t.__origin__, List)):
            return list(flatten([flatten_type(sub_t) for sub_t in t.__args__]))
    return [t]


def overlaps(a: list, b: list):
    """
    Returns true if there are any elements in common between a and b
    """
    return len(set(a).intersection(set(b))) > 0


def pydantic_class_refatt_map():
    """
    Builds a map of class names to their field names that are referable types.
    As in, types with an identifier that can be referred to elsewhere,
    collapsed to that identifier and dereferenced.
    Returns a map like:
    {"Allele": ["location"], ...}
    """
    # Things defined here that are classes that inherit from BaseModel
    this_module = sys.modules[__name__]
    global_map = globals()
    model_classes = list(filter(
        lambda c: (inspect.isclass(c)
                   and issubclass(c, BaseModel)
                   and inspect.getmodule(c) == this_module),
        [gl_name_value[1] for gl_name_value in global_map.items()]
    ))
    # Types directly reffable
    reffable_classes = list(filter(
        lambda c: ('id' in c.model_fields
                   and c.is_ga4gh_identifiable()),
        model_classes
    ))
    # Types reffable because they are a union of reffable types
    union_reffable_classes = []
    for model_class in model_classes:
        if issubclass(model_class, RootModel):
            flattened_type_annotation = flatten_type(model_class.model_fields["root"].annotation)
            if overlaps(reffable_classes, flattened_type_annotation):
                union_reffable_classes.append(model_class)
    reffable_fields = {}
    # Find any field whose type is a subclass of a reffable type,
    # or which is a typing.List that includes a reffable type.
    for model_class in model_classes:
        fields = model_class.model_fields
        class_reffable_fields = []
        for fieldname, field in fields.items():
            if fieldname == 'root':
                continue
            field_type = field.annotation  # a typing or normal annotation like str
            # types can be raw class type annotation (int, str, dict, etc)
            # typing.Literal, typing.Union, typing.Optional
            # Use flatten_type to simplify these
            if any([rc in flatten_type(field_type)
                    for rc in (reffable_classes
                               + union_reffable_classes)]):
                class_reffable_fields.append(fieldname)
        if len(class_reffable_fields) > 0:
            reffable_fields[model_class.__name__] = class_reffable_fields
    class_keys = {}
    for model_class in model_classes:
        keys = getattr_in(model_class, ['ga4gh', 'keys'])
        if keys and len(keys) > 0:
            class_keys[model_class.__name__] = keys
    return (reffable_classes,
            union_reffable_classes,
            reffable_fields,
            class_keys)


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
    SEQ_TERMINUS = "SequenceTerminus"
    DERIVATIVE_SEQ = "DerivativeSequence"
    CN_COUNT = "CopyNumberCount"
    CN_CHANGE = "CopyNumberChange"


class ResidueAlphabet(str, Enum):
    """The interpretation of the character codes referred to by the refget accession,
    where "aa" specifies an amino acid character set, and "na" specifies a nucleic acid
    character set.
    """

    AA = 'aa'
    NA = 'na'


class CopyChange(str, Enum):
    """Define constraints for copy change"""

    EFO_0030069 = 'efo:0030069'
    EFO_0020073 = 'efo:0020073'
    EFO_0030068 = 'efo:0030068'
    EFO_0030067 = 'efo:0030067'
    EFO_0030064 = 'efo:0030064'
    EFO_0030070 = 'efo:0030070'
    EFO_0030071 = 'efo:0030071'
    EFO_0030072 = 'efo:0030072'


def _recurse_ga4gh_serialize(obj):
    if isinstance(obj, _Ga4ghIdentifiableObject):
        return obj.get_or_create_digest()
    elif isinstance(obj, _ValueObject):
        return obj.ga4gh_serialize()
    elif isinstance(obj, RootModel):
        return _recurse_ga4gh_serialize(obj.model_dump())
    elif isinstance(obj, str):
        return obj
    elif isinstance(obj, list):
        return [_recurse_ga4gh_serialize(x) for x in obj]
    else:
        return obj


class _ValueObject(DomainEntity, ABC):
    """A contextual value whose equality is based on value, not identity.
    See https://en.wikipedia.org/wiki/Value_object for more on Value Objects.
    """

    def __hash__(self):
        return encode_canonical_json(self.ga4gh_serialize()).decode("utf-8").__hash__()

    def ga4gh_serialize(self) -> Dict:
        out = OrderedDict()
        for k in self.ga4gh.keys:
            v = getattr(self, k)
            out[k] = _recurse_ga4gh_serialize(v)
        return out

    class ga4gh:
        keys: List[str]

    @staticmethod
    def is_ga4gh_identifiable():
        return False


class _Ga4ghIdentifiableObject(_ValueObject, ABC):
    """A contextual value object for which a GA4GH computed identifier can be created.
    All GA4GH Identifiable Objects may have computed digests from the VRS Computed
    Identifier algorithm.
    """

    type: str
    digest: Optional[Annotated[str, StringConstraints(pattern=r'^[0-9A-Za-z_\-]{32}$')]] = Field(
        None,
        description='A sha512t24u digest created using the VRS Computed Identifier algorithm.',
    )

    def __lt__(self, other):
        return self.get_or_create_digest() < other.get_or_create_digest()

    @staticmethod
    def is_ga4gh_identifiable():
        return True

    def has_valid_ga4gh_id(self):
        return self.id and GA4GH_IR_REGEXP.match(self.id) is not None

    def compute_digest(self, store=True, as_version: PrevVrsVersion | None = None) -> str:
        """A sha512t24u digest created using the VRS Computed Identifier algorithm.

        Stores the digest in the object if ``store`` is ``True``.

        If ``as_version`` is provided, other parameters are ignored and a digest is
        returned following the conventions of the VRS version indicated by ``as_version_``.
        """
        if as_version is None:
            digest = sha512t24u(encode_canonical_json(self.ga4gh_serialize()))
            if store:
                self.digest = digest
        else:
            try:
                digest = sha512t24u(self.ga4gh_serialize_as_version(as_version).encode("utf-8"))
            except AttributeError:
                raise AttributeError('This class does not support prior version identifiers.')
        return digest

    def get_or_create_ga4gh_identifier(self, in_place='default', recompute=False, as_version=None) -> str:
        """Sets and returns a GA4GH Computed Identifier for the object.
        Overwrites the existing identifier if overwrite is True.

        This function has three options for in_place editing of vro.id:
        - 'default': the standard identifier update behavior for GA4GH
            identifiable objects, this mode will update the vro.id
            field if the field is empty
        - 'always': this will update the vro.id field any time the
            identifier is computed
        - 'never': the vro.id field will not be edited in-place,
            even when empty

        Digests will be recalculated even if present if recompute is True.

        If ``as_version`` is provided, other parameters are ignored and an identifier is
        returned following the conventions of the VRS version indicated by
        ``as_version_``.
        """
        if as_version is not None:
            return self.compute_ga4gh_identifier(as_version=as_version)

        if in_place == 'default':
            if self.id is None:
                self.id = self.compute_ga4gh_identifier(recompute)
        elif in_place == 'always':
            self.id = self.compute_ga4gh_identifier(recompute)
        elif in_place == 'never':
            return self.compute_ga4gh_identifier(recompute)
        else:
            raise ValueError("Expected 'in_place' to be one of 'default', 'always', or 'never'")

        if self.has_valid_ga4gh_id():
            return self.id
        else:
            return self.compute_ga4gh_identifier(recompute)

    def compute_ga4gh_identifier(self, recompute=False, as_version=None):
        """Returns a GA4GH Computed Identifier.

        If ``as_version`` is provided, other parameters are ignored and a computed
        identifier is returned following the conventions of the VRS version indicated by
        ``as_version_``.
        """
        if as_version is None:
            self.get_or_create_digest(recompute)
            return f'{CURIE_NAMESPACE}{CURIE_SEP}{self.ga4gh.prefix}{GA4GH_PREFIX_SEP}{self.digest}'
        else:
            digest = self.compute_digest(as_version=as_version)
            return f'{CURIE_NAMESPACE}{CURIE_SEP}{self.ga4gh.priorPrefix[as_version]}{GA4GH_PREFIX_SEP}{digest}'

    def get_or_create_digest(self, recompute=False) -> str:
        """Sets and returns a sha512t24u digest of the GA4GH Identifiable Object, or creates
        the digest if it does not exist."""
        if self.digest is None or recompute:
            return self.compute_digest()
        return self.digest

    class ga4gh(_ValueObject.ga4gh):
        prefix: str


#########################################
# vrs numerics, comparators, and ranges
#########################################

class Range(RootModel):
    """An inclusive range of values bounded by one or more integers."""

    root: List[Optional[int]] = Field(
        ...,
        json_schema_extra={
            'description': 'An inclusive range of values bounded by one or more integers.'
        },
        max_length=2,
        min_length=2,
    )


class Residue(RootModel):
    """A character representing a specific residue (i.e., molecular species) or
    groupings of these ("ambiguity codes"), using `one-letter IUPAC abbreviations
    <https://en.wikipedia.org/wiki/International_Union_of_Pure_and_Applied_Chemistry#Amino_acid_and_nucleotide_base_codes>`_
    for nucleic acids and amino acids.
    """

    root: Annotated[str, StringConstraints(pattern=r'[A-Z*\-]')] = Field(
        ...,
        json_schema_extra={
            'description': 'A character representing a specific residue (i.e., molecular species) or groupings of these ("ambiguity codes"), using [one-letter IUPAC abbreviations](https://en.wikipedia.org/wiki/International_Union_of_Pure_and_Applied_Chemistry#Amino_acid_and_nucleotide_base_codes) for nucleic acids and amino acids.'
        },
    )


class SequenceString(RootModel):
    """A character string of `Residues` that represents a biological sequence using the
    conventional sequence order (5'-to-3' for nucleic acid sequences, and
    amino-to-carboxyl for amino acid sequences). IUPAC ambiguity codes are permitted in
    Sequence Strings.
    """

    root: Annotated[str, StringConstraints(pattern=r'^[A-Z*\-]*$')] = Field(
        ...,
        json_schema_extra={
            'description': "A character string of Residues that represents a biological sequence using the conventional sequence order (5'-to-3' for nucleic acid sequences, and amino-to-carboxyl for amino acid sequences). IUPAC ambiguity codes are permitted in Sequence Strings."
        },
    )


#########################################
# vrs sequence expression
#########################################


class LengthExpression(_ValueObject):
    """A sequence expressed only by its length."""

    type: Literal["LengthExpression"] = Field(
        VrsType.LEN_EXPR.value, description=f'MUST be "{VrsType.LEN_EXPR.value}"'
    )
    length: Optional[Union[Range, int]] = None

    class ga4gh(_ValueObject.ga4gh):
        keys = [
            'length',
            'type'
        ]


class ReferenceLengthExpression(_ValueObject):
    """An expression of a length of a sequence from a repeating reference."""

    type: Literal["ReferenceLengthExpression"] = Field(
        VrsType.REF_LEN_EXPR.value, description=f'MUST be "{VrsType.REF_LEN_EXPR.value}"'
    )
    length: Union[Range, int] = Field(
        ..., description='The number of residues of the expressed sequence.'
    )
    sequence: Optional[SequenceString] = Field(
        None, description='the `Sequence` encoded by the Reference Length Expression.'
    )
    repeatSubunitLength: int = Field(
        ..., description='The number of residues of the repeat subunit.'
    )

    class ga4gh(_ValueObject.ga4gh):
        keys = [
            'length',
            'repeatSubunitLength',
            'type'
        ]


class LiteralSequenceExpression(_ValueObject):
    """An explicit expression of a Sequence."""

    type: Literal["LiteralSequenceExpression"] = Field(
        VrsType.LIT_SEQ_EXPR.value, description=f'MUST be "{VrsType.LIT_SEQ_EXPR.value}"'
    )
    sequence: SequenceString = Field(..., description='the literal sequence')

    class ga4gh(_ValueObject.ga4gh):
        keys = [
            'sequence',
            'type'
        ]


#########################################
# vrs location
#########################################


class SequenceReference(_ValueObject):
    """A sequence of nucleic or amino acid character codes."""

    model_config = ConfigDict(use_enum_values=True)

    type: Literal["SequenceReference"] = Field(VrsType.SEQ_REF.value, description=f'MUST be "{VrsType.SEQ_REF.value}"')
    refgetAccession: Annotated[str, StringConstraints(pattern=r'^SQ.[0-9A-Za-z_\-]{32}$')] = Field(
        ...,
        description='A `GA4GH RefGet <http://samtools.github.io/hts-specs/refget.html>` identifier for the referenced sequence, using the sha512t24u digest.',
    )
    residueAlphabet: Optional[ResidueAlphabet] = Field(None, description="The interpretation of the character codes referred to by the refget accession, where 'aa' specifies an amino acid character set, and 'na' specifies a nucleic acid character set.")
    circular: Optional[bool] = Field(None, description="A boolean indicating whether the molecule represented by the sequence is circular (true) or linear (false).")

    class ga4gh(_ValueObject.ga4gh):
        keys = [
            'refgetAccession',
            'type'
        ]


class SequenceLocation(_Ga4ghIdentifiableObject):
    """A `Location` defined by an interval on a referenced `Sequence`."""

    type: Literal["SequenceLocation"] = Field(VrsType.SEQ_LOC.value, description=f'MUST be "{VrsType.SEQ_LOC.value}"')
    sequenceReference: Optional[Union[IRI, SequenceReference]] = Field(
        None, description='A reference to a `Sequence` on which the location is defined.'
    )
    start: Optional[Union[Range, int]] = Field(
        None,
        description='The start coordinate or range of the SequenceLocation. The minimum value of this coordinate or range is 0. MUST represent a coordinate or range less than the value of `end`.',
    )
    end: Optional[Union[Range, int]] = Field(
        None,
        description='The end coordinate or range of the SequenceLocation. The minimum value of this coordinate or range is 0. MUST represent a coordinate or range greater than the value of `start`.',

    )
    sequence: Optional[SequenceString] = Field(None, description="The literal sequence encoded by the `sequenceReference` at these coordinates.")

    def ga4gh_serialize_as_version(self, as_version: PrevVrsVersion):
        """This method will return a serialized string following the conventions for
        SequenceLocation serialization as defined in the VRS version specified by
        ``as_version``.

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
                    raise ValueError(f'{value} is not int or list.')
                out.append(result)
            return f'{{"interval":{{"end":{out[1]},"start":{out[0]},"type":"SequenceInterval"}},"sequence_id":"{self.sequenceReference.refgetAccession.split(".")[1]}","type":"SequenceLocation"}}'

    def get_refget_accession(self):
        if isinstance(self.sequenceReference, SequenceReference):
            return self.sequenceReference.refgetAccession
        elif isinstance(self.sequenceReference, IRI):
            return self.sequenceReference.root
        else:
            return None

    class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
        prefix = 'SL'
        priorPrefix = {PrevVrsVersion.V1_3.value: 'VSL'}
        keys = [
            'end',
            'sequenceReference',
            'start',
            'type'
        ]

#########################################
# base variation
#########################################


class _VariationBase(_Ga4ghIdentifiableObject, ABC):
    """Base class for variation"""

    expressions: Optional[List[Expression]] = None

#########################################
# vrs molecular variation
#########################################


class Allele(_VariationBase):
    """The state of a molecule at a `Location`."""

    type: Literal["Allele"] = Field(VrsType.ALLELE.value, description=f'MUST be "{VrsType.ALLELE.value}"')
    location: Union[IRI, SequenceLocation] = Field(
        ..., description='The location of the Allele'
    )
    state: Union[LiteralSequenceExpression, ReferenceLengthExpression, LengthExpression] = Field(
        ..., description='An expression of the sequence state'
    )

    def ga4gh_serialize_as_version(self, as_version: PrevVrsVersion):
        """This method will return a serialized string following the conventions for
        Allele serialization as defined in the VRS version specified by 'as_version`.

        :raises ValueError: If ``state`` is not a ``LiteralSequenceExpression`` or
            ``ReferenceLengthExpression``; ``state.sequence`` is null.
        """
        location_digest = self.location.compute_digest(as_version=as_version)

        if not isinstance(self.state, (LiteralSequenceExpression, ReferenceLengthExpression)):
            err_msg = "Only `LiteralSequenceExpression` and `ReferenceLengthExpression` are supported for previous versions of VRS"
            raise ValueError(err_msg)

        sequence = get_pydantic_root(self.state.sequence)

        if sequence is None:
            raise ValueError('State sequence attribute must be defined.')

        if as_version == PrevVrsVersion.V1_3:
            return f'{{"location":"{location_digest}","state":{{"sequence":"{sequence}","type":"LiteralSequenceExpression"}},"type":"Allele"}}'


    class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
        prefix = 'VA'
        priorPrefix = {PrevVrsVersion.V1_3.value: 'VA'}
        keys = [
            'location',
            'state',
            'type'
        ]


class CisPhasedBlock(_VariationBase):
    """An ordered set of co-occurring `Variation` on the same molecule."""

    type: Literal["CisPhasedBlock"] = Field(VrsType.CIS_PHASED_BLOCK.value, description=f'MUST be "{VrsType.CIS_PHASED_BLOCK.value}"')
    members: List[Union[Allele, IRI]] = Field(
        ...,
        description='A list of `Alleles` that are found in-cis on a shared molecule.',
        min_length=2,
    )
    sequenceReference: Optional[SequenceReference] = Field(None, description="An optional Sequence Reference on which all of the in-cis Alleles are found. When defined, this may be used to implicitly define the `sequenceReference` attribute for each of the CisPhasedBlock member Alleles.")

    def ga4gh_serialize(self) -> Dict:
        out = _ValueObject.ga4gh_serialize(self)
        out["members"] = sorted(out["members"])
        return out

    class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
        prefix = 'CPB'
        keys = [
            'members',
            'type'
        ]


#########################################
# vrs structural variation (under active discussion)
#########################################


class Adjacency(_VariationBase):
    """The `Adjacency` class can represent either the termination of a sequence or the
    adjoining of the end of a sequence with the beginning of an adjacent sequence,
    potentially with an intervening linker sequence.
    """

    type: Literal["Adjacency"] = Field(VrsType.ADJACENCY.value, description=f'MUST be "{VrsType.ADJACENCY.value}"')
    adjoinedSequences: List[Union[IRI, SequenceLocation]] = Field(
        ...,
        description="The terminal sequence or pair of adjoined sequences that defines in the adjacency.",
        min_length=2,
        max_length=2,
    )
    linker: Optional[Union[LiteralSequenceExpression, ReferenceLengthExpression, LengthExpression]] = Field(
        None,
        description="The sequence found between adjoined sequences."
    )
    homology: Optional[bool] = Field(None, description="A flag indicating if coordinate ambiguity in the adjoined sequences is from sequence homology (true) or other uncertainty (false).")

    class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
        prefix = 'AJ'
        keys = [
            'adjoinedSequences',
            'linker',
            'type'
        ]


class SequenceTerminus(_VariationBase):
    """The `SequenceTerminus` data class provides a structure for describing the end
    (terminus) of a sequence. Structurally similar to Adjacency but the linker sequence
    is not allowed and it removes the unnecessary array structure.
    """

    type: Literal["SequenceTerminus"] = Field(VrsType.SEQ_TERMINUS.value, description=f'MUST be "{VrsType.SEQ_TERMINUS.value}"')
    location: Union[IRI, SequenceLocation] = Field(..., description="The location of the terminus.")

    class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
        prefix = "SQX"
        keys = [
            "location",
            "type"
        ]


class DerivativeSequence(_VariationBase):
    """The "Derivative Sequence" data class is a structure for describing a derivate
    sequence composed from multiple sequence adjacencies.
    """

    type: Literal["DerivativeSequence"] = Field(VrsType.DERIVATIVE_SEQ.value, description=f'MUST be "{VrsType.DERIVATIVE_SEQ.value}"')
    components: List[Union[IRI, Adjacency, Allele, SequenceTerminus, CisPhasedBlock]] = Field(
        ...,
        description="The sequence components that make up the derivative sequence.",
        min_length=2
    )

    class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
        prefix = "DSQ"
        keys = [
            "components",
            "type"
        ]


#########################################
# vrs systemic variation
#########################################


class _CopyNumber(_VariationBase, ABC):
    """A measure of the copies of a `Location` within a system (e.g. genome, cell, etc.)"""

    location: Union[IRI, SequenceLocation] = Field(
        ...,
        description='A location for which the number of systemic copies is described.',
    )


class CopyNumberCount(_CopyNumber):
    """The absolute count of discrete copies of a `Location` or `Gene`, within a system
    (e.g. genome, cell, etc.).
    """

    type: Literal["CopyNumberCount"] = Field(VrsType.CN_COUNT.value, description=f'MUST be "{VrsType.CN_COUNT.value}"')
    copies: Union[Range, int] = Field(
        ..., description='The integral number of copies of the subject in a system'
    )

    class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
        prefix = 'CN'
        keys = [
            'copies',
            'location',
            'type'
        ]


class CopyNumberChange(_CopyNumber):
    """An assessment of the copy number of a `Location` or a `Gene` within a system
    (e.g. genome, cell, etc.) relative to a baseline ploidy.
    """

    model_config = ConfigDict(use_enum_values=True)

    type: Literal["CopyNumberChange"] = Field(VrsType.CN_CHANGE.value, description=f'MUST be "{VrsType.CN_CHANGE.value}"')
    copyChange: CopyChange = Field(
        ...,
        description='MUST be one of "efo:0030069" (complete genomic loss), "efo:0020073" (high-level loss), "efo:0030068" (low-level loss), "efo:0030067" (loss), "efo:0030064" (regional base ploidy), "efo:0030070" (gain), "efo:0030071" (low-level gain), "efo:0030072" (high-level gain).',
    )

    class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
        prefix = 'CX'
        keys = [
            'copyChange',
            'location',
            'type'
        ]


#########################################
# vrs kinds of variation, expression, and location
#########################################


class MolecularVariation(RootModel):
    """A `variation` on a contiguous molecule."""

    root: Union[Allele, CisPhasedBlock, Adjacency, SequenceTerminus, DerivativeSequence] = Field(
        ...,
        json_schema_extra={
            'description': 'A `variation` on a contiguous molecule.'
        },
        discriminator='type'
    )

class SequenceExpression(RootModel):
    """An expression describing a `Sequence`."""

    root: Union[LiteralSequenceExpression, ReferenceLengthExpression, LengthExpression] = Field(
        ...,
        json_schema_extra={'description': 'An expression describing a `Sequence`.'},
        discriminator='type'
    )


class Location(RootModel):
    """A contiguous segment of a biological sequence."""

    root: SequenceLocation = Field(
        ...,
        json_schema_extra={
            'description': 'A contiguous segment of a biological sequence.'
        },
        discriminator='type',
    )


class Variation(RootModel):
    """A representation of the state of one or more biomolecules."""

    root: Union[Allele, CisPhasedBlock, Adjacency, SequenceTerminus, DerivativeSequence, CopyNumberChange, CopyNumberCount] = Field(
        ...,
        json_schema_extra={
            'description': 'A representation of the state of one or more biomolecules.'
        },
        discriminator='type',
    )


class SystemicVariation(RootModel):
    """A Variation of multiple molecules in the context of a system, e.g. a genome,
    sample, or homologous chromosomes.
    """

    root: Union[CopyNumberChange, CopyNumberCount] = Field(
        ...,
        json_schema_extra={
            'description': 'A Variation of multiple molecules in the context of a system, e.g. a genome, sample, or homologous chromosomes.'
        },
        discriminator='type',
    )


# At end so classes exist
(
    reffable_classes,
    union_reffable_classes,
    class_refatt_map,
    class_keys
) = pydantic_class_refatt_map()
