"""
**This module should not be imported directly.**

Instead, users should use one of the following:

  * `from ga4gh.vrs import models`, and refer to models with the
    abbreviated name, e.g., `models.Allele` (recommended)

  * `import ga4gh.vrs`, and refer to models using the fully-qualified
    module name, e.g., `ga4gh.vrs.models.Allele`

New pydantic-based version

Pydantic classes bootstrapped with:
sed -i.bkp 's/$defs/definitions/g' merged.json
V1 pydantic: datamodel-codegen --input submodules/vrs/schema/merged.json --input-file-type jsonschema --output models_merged2.py
V2 pydantic: datamodel-codegen --input submodules/vrs/schema/merged.json --input-file-type jsonschema --output models.py --output-model-type pydantic_v2.BaseModel --allow-extra-fields
"""

from typing import List, Literal, Optional, Union, Dict, Set
from collections import OrderedDict
from enum import Enum
import inspect
import sys
import typing
from ga4gh.core import sha512t24u, GA4GH_PREFIX_SEP, CURIE_SEP, CURIE_NAMESPACE, GA4GH_IR_REGEXP

from pydantic import BaseModel, ConfigDict, Field, RootModel, constr, model_serializer

from ga4gh.core._internal.pydantic import (
    is_ga4gh_identifiable,
    getattr_in
)
from ga4gh.core._internal.models import IRI, _Entity


def flatten(vals):
    """
    Flattens vals recursively, lazily using yield
    """
    def is_coll(thing):
        """
        Return True if the thing looks like a collection.
        This is not exhaustive, do not use in general.
        """
        # return hasattr(thing, '__iter__') and not isinstance(thing, str) and not inspect.isclass(thing)
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
        if t.__origin__ == typing.Literal:
            return list(t.__args__)
        elif (t.__origin__ == typing.Union
              or issubclass(t.__origin__, typing.List)):
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
                   and is_ga4gh_identifiable(c)),
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


class Syntax(Enum):
    """Define constraints for syntax"""

    HGVS_C = "hgvs.c"
    HGVS_P = "hgvs.p"
    HGVS_G = "hgvs.g"
    HGVS_M = "hgvs.m"
    HGVS_N = "hgvs.n"
    HGVS_R = "hgvs.r"
    ISCN = "iscn"
    GNOMAD = "gnomad"
    SPDI = "spdi"


class ResidueAlphabet(Enum):
    AA = 'aa'
    NA = 'na'


class CopyChange(Enum):
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
        if obj.startswith(f'{CURIE_NAMESPACE}{CURIE_SEP}'):
            return obj.split(GA4GH_PREFIX_SEP)[-1]
        return obj
    # This is code that may be used to handle unordered lists (sets) as defined
    # by VRS. This functionality will be needed for Haplotype if we return to ordered: false.
    # No attempt is made to order arrays of strings containing objects. That will need
    # to be addressed before implementing the non-identifiable GenotypeMember class.
    elif isinstance(obj, set):
        out = [_recurse_ga4gh_serialize(x) for x in list(obj)]
        if all(isinstance(x, str) for x in out):
            return sorted(out)
        else:
            return out
    elif isinstance(obj, list):
        return [_recurse_ga4gh_serialize(x) for x in obj]
    else:
        return obj


class _ValueObject(_Entity):
    """A contextual value whose equality is based on value, not identity.
    See https://en.wikipedia.org/wiki/Value_object for more on Value Objects.
    """

    def __hash__(self):
        return self.model_dump_json().__hash__()

    @model_serializer(when_used='json')
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


class _Ga4ghIdentifiableObject(_ValueObject):
    """A contextual value object for which a GA4GH computed identifier can be created.
    All GA4GH Identifiable Objects may have computed digests from the VRS Computed
    Identifier algorithm."""

    type: str

    digest: Optional[constr(pattern=r'^[0-9A-Za-z_\-]{32}$')] = Field(
        None,
        description='A sha512t24u digest created using the VRS Computed Identifier algorithm.',
    )

    @staticmethod
    def is_ga4gh_identifiable():
        return True

    def has_valid_ga4gh_id(self):
        return self.id and GA4GH_IR_REGEXP.match(self.id) is not None

    def has_valid_digest(self):
        return bool(self.digest)  # Pydantic constraint ensures digest field value is valid

    def compute_digest(self, store=True) -> str:
        """A sha512t24u digest created using the VRS Computed Identifier algorithm.
        Stores the digest in the object if store is True.
        """
        digest = sha512t24u(self.model_dump_json().encode("utf-8"))
        if store:
            self.digest = digest
        return digest

    def get_or_create_ga4gh_identifier(self, overwrite=False) -> str:
        """Sets and returns a GA4GH Computed Identifier for the object.
        Overwrites the existing identifier if overwrite is True."""
        if self.id is None or overwrite:
            self.get_or_create_digest()
            self.id = f'{CURIE_NAMESPACE}{CURIE_SEP}{self.ga4gh.prefix}{GA4GH_PREFIX_SEP}{self.digest}'
        return self.id

    def get_or_create_digest(self, overwrite=False) -> str:
        """Sets and returns a sha512t24u digest of the GA4GH Identifiable Object, or creates
        the digest if it does not exist."""
        if self.digest is None or overwrite:
            return self.compute_digest()
        return self.digest

    class ga4gh(_ValueObject.ga4gh):
        prefix: str


class Expression(BaseModel):
    """Representation of a variation by a specified nomenclature or syntax for a
    Variation object. Common examples of expressions for the description of molecular
    variation include the HGVS and ISCN nomenclatures.
    """
    model_config = ConfigDict(
        use_enum_values=True
    )

    syntax: Syntax
    value: str
    syntax_version: Optional[str] = None


class Range(RootModel):
    root: List[Optional[int]] = Field(
        ...,
        json_schema_extra={
            'description': 'An inclusive range of values bounded by one or more integers.'
        },
        max_length=2,
        min_length=2,
    )


class Residue(RootModel):
    root: constr(pattern=r'[A-Z*\-]') = Field(
        ...,
        json_schema_extra={
            'description': 'A character representing a specific residue (i.e., molecular species) or groupings of these ("ambiguity codes"), using [one-letter IUPAC abbreviations](https://en.wikipedia.org/wiki/International_Union_of_Pure_and_Applied_Chemistry#Amino_acid_and_nucleotide_base_codes) for nucleic acids and amino acids.'
        },
    )


class SequenceString(RootModel):
    root: constr(pattern=r'^[A-Z*\-]*$') = Field(
        ...,
        json_schema_extra={
            'description': 'A character string of Residues that represents a biological sequence using the conventional sequence order (5’-to-3’ for nucleic acid sequences, and amino-to-carboxyl for amino acid sequences). IUPAC ambiguity codes are permitted in Sequence Strings.'
        },
    )


class SequenceReference(_ValueObject):
    model_config = ConfigDict(
        use_enum_values=True
    )

    type: Literal['SequenceReference'] = Field('SequenceReference', description='MUST be "SequenceReference"')
    refgetAccession: constr(pattern=r'^SQ.[0-9A-Za-z_\-]{32}$') = Field(
        ...,
        description='A `GA4GH RefGet <http://samtools.github.io/hts-specs/refget.html>` identifier for the referenced sequence, using the sha512t24u digest.',
    )
    residueAlphabet: Optional[ResidueAlphabet] = None

    class ga4gh(_ValueObject.ga4gh):
        keys = [
            'refgetAccession',
            'type'
        ]


class LengthExpression(_ValueObject):
    """An expression of a DNA, RNA, or protein polymer of known length but unspecified sequence."""

    type: Literal['ReferenceLengthExpression'] = Field(
        'ReferenceLengthExpression', description='MUST be "ReferenceLengthExpression"'
    )
    length: Union[Range, int] = Field(
        ..., description='The number of residues of the expressed sequence.'
    )

    class ga4gh(_ValueObject.ga4gh):
        keys = [
            'length',
            'type'
        ]


class ReferenceLengthExpression(_ValueObject):
    """An expression sequence derived from a reference."""

    type: Literal['ReferenceLengthExpression'] = Field(
        'ReferenceLengthExpression', description='MUST be "ReferenceLengthExpression"'
    )
    length: Union[Range, int] = Field(
        ..., description='The number of residues of the expressed sequence.'
    )
    sequence: Optional[SequenceString] = Field(
        None, description='the Sequence encoded by the Reference Length Expression.'
    )
    repeatSubunitLength: int = Field(
        None, description='The number of residues of the repeat subunit.'
    )

    class ga4gh(_ValueObject.ga4gh):
        keys = [
            'length',
            'repeatSubunitLength',
            'type'
        ]


class LiteralSequenceExpression(_ValueObject):
    """An explicit expression of a Sequence."""

    type: Literal['LiteralSequenceExpression'] = Field(
        'LiteralSequenceExpression', description='MUST be "LiteralSequenceExpression"'
    )
    sequence: SequenceString = Field(..., description='the literal sequence')

    class ga4gh(_ValueObject.ga4gh):
        keys = [
            'sequence',
            'type'
        ]


class SequenceLocation(_Ga4ghIdentifiableObject):
    """A `Location` defined by an interval on a referenced `Sequence`."""

    type: Literal['SequenceLocation'] = Field('SequenceLocation', description='MUST be "SequenceLocation"')
    sequenceReference: Optional[Union[IRI, SequenceReference]] = Field(
        None, description='A SequenceReference.'
    )
    start: Union[Range, int] = Field(
        ...,
        description='The start coordinate or range of the SequenceLocation. The minimum value of this coordinate or range is 0. MUST represent a coordinate or range less than the value of `end`.',
    )
    end: Union[Range, int] = Field(
        ...,
        description='The end coordinate or range of the SequenceLocation. The minimum value of this coordinate or range is 0. MUST represent a coordinate or range greater than the value of `start`.',
    )

    class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
        prefix = 'SL'
        keys = [
            'end',
            'sequenceReference',
            'start',
            'type'
        ]


class _VariationBase(_Ga4ghIdentifiableObject):
    """Base class for variation"""

    expressions: Optional[List[Expression]] = None


class Allele(_VariationBase):

    type: Literal['Allele'] = Field('Allele', description='MUST be "Allele"')
    location: Union[IRI, SequenceLocation] = Field(
        ..., description='The location of the Allele'
    )
    state: Union[LiteralSequenceExpression, ReferenceLengthExpression] = Field(
        ..., description='An expression of the sequence state'
    )

    class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
        prefix = 'VA'
        keys = [
            'location',
            'state',
            'type'
        ]


class Haplotype(_VariationBase):
    """A set of non-overlapping Allele members that co-occur on the same molecule."""

    type: Literal['Haplotype'] = Field('Haplotype', description='MUST be "Haplotype"')
    # TODO members temporarily typed as Set instead of List
    members: Set[Union[Allele, IRI]] = Field(
        ...,
        description='A list of Alleles (or IRI references to `Alleles`) that comprise a Haplotype. Since each `Haplotype` member MUST be an `Allele`, and all members MUST share a common `SequenceReference`, implementations MAY use a compact representation of Haplotype that omits type and `SequenceReference` information in individual Haplotype members. Implementations MUST transform compact `Allele` representations into an `Allele` when computing GA4GH identifiers.',
        min_length=2,
    )

    class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
        prefix = 'HT'
        keys = [
            'members',
            'type'
        ]


class _CopyNumber(_VariationBase):
    """A measure of the copies of a `Location` within a system (e.g. genome, cell, etc.)"""

    location: Union[IRI, SequenceLocation] = Field(
        ...,
        description='A location for which the number of systemic copies is described.',
    )


class CopyNumberCount(_CopyNumber):
    """The absolute count of discrete copies of a `Location` or `Gene`, within a system
    (e.g. genome, cell, etc.).
    """

    type: Literal['CopyNumberCount'] = Field('CopyNumberCount', description='MUST be "CopyNumberCount"')
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
    model_config = ConfigDict(
        use_enum_values=True
    )

    type: Literal['CopyNumberChange'] = Field('CopyNumberChange', description='MUST be "CopyNumberChange"')
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


# class GenotypeMember(_ValueObject):
#     """A class for expressing the count of a specific `MolecularVariation` present
#     in-trans at a genomic locus represented by a `Genotype`.
#     """
#
#     type: Literal['GenotypeMember'] = Field('GenotypeMember', description='MUST be "GenotypeMember".')
#     count: Union[Range, int] = Field(
#         ..., description='The number of copies of the `variation` at a Genotype locus.'
#     )
#     variation: Union[Allele, Haplotype] = Field(
#         ..., description='A MolecularVariation at a Genotype locus.'
#     )
#
#     class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
#         keys = [
#             'type',
#             'count',
#             'variation'
#         ]


class MolecularVariation(RootModel):
    """A variation on a contiguous molecule."""

    root: Union[Allele, Haplotype] = Field(
        ...,
        json_schema_extra={
            'description': 'A variation on a contiguous molecule.'
        },
        discriminator='type'
    )


# class Genotype(_VariationBase):
#     """A quantified set of _in-trans_ `MolecularVariation` at a genomic locus."""
#
#     type: Literal['Genotype'] = Field(
#         'Genotype',
#         description='MUST be "Genotype"'
#     )
#     # TODO members temporarily typed as List instead of Set + validate unique items
#     members: List[GenotypeMember] = Field(
#         ...,
#         description='Each GenotypeMember in `members` describes a MolecularVariation and the count of that variation at the locus.',
#         min_length=1,
#     )
#     count: Union[Range, int] = Field(
#         ...,
#         description='The total number of copies of all MolecularVariation at this locus, MUST be greater than or equal to the sum of GenotypeMember copy counts. If greater than the total counts, this implies additional MolecularVariation that are expected to exist but are not explicitly indicated.',
#     )
#
#     class ga4gh(_Ga4ghIdentifiableObject.ga4gh):
#         prefix = 'GT'
#         keys = [
#             'count',
#             'members',
#             'type'
#         ]


class SequenceExpression(RootModel):
    root: Union[LiteralSequenceExpression, ReferenceLengthExpression] = Field(
        ...,
        json_schema_extra={'description': 'An expression describing a Sequence.'},
        discriminator='type'
    )


class Location(RootModel):
    root: SequenceLocation = Field(
        ...,
        json_schema_extra={
            'description': 'A contiguous segment of a biological sequence.'
        },
        discriminator='type',
    )


class Variation(RootModel):
    root: Union[Allele, CopyNumberChange, CopyNumberCount, Haplotype] = Field(
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
