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

from typing import Any, Dict, List, Literal, Optional, Union
from enum import Enum
import inspect
import logging
import sys
import typing

from pydantic import BaseModel, ConfigDict, Field, RootModel, constr, model_validator

from ga4gh.core._internal.pydantic import (
    is_identifiable,
    getattr_in
)

_logger = logging.getLogger(__name__)


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
                   and is_identifiable(c)),
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


class Relation(Enum):
    """A mapping relation between concepts as defined by the Simple Knowledge
    Organization System (SKOS).
    """

    CLOSE_MATCH = 'closeMatch'
    EXACT_MATCH = 'exactMatch'
    BROAD_MATCH = 'broadMatch'
    NARROW_MATCH = 'narrowMatch'
    RELATED_MATCH = 'relatedMatch'


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


class Code(RootModel):
    """Indicates that the value is taken from a set of controlled strings defined
    elsewhere. Technically, a code is restricted to a string which has at least one
    character and no leading or  trailing whitespace, and where there is no whitespace
    other than single spaces in the contents."""

    root: constr(pattern=r'\S+( \S+)*') = Field(
        ...,
        description='Indicates that the value is taken from a set of controlled strings defined elsewhere. Technically, a code is restricted to a string which has at least one character and no leading or  trailing whitespace, and where there is no whitespace other than single spaces in the contents.',
        example='ENSG00000139618',
    )


class IRI(RootModel):
    root: str = Field(
        ...,
        description='An IRI Reference (either an IRI or a relative-reference), according to `RFC3986 section 4.1  <https://datatracker.ietf.org/doc/html/rfc3986#section-4.1>` and `RFC3987 section 2.1 <https://datatracker.ietf.org/doc/html/rfc3987#section-2.1>`. MAY be a JSON Pointer as an IRI fragment, as  described by `RFC6901 section 6 <https://datatracker.ietf.org/doc/html/rfc6901#section-6>`.',
    )


class Extension(BaseModel):
    """The Extension class provides VODs with a means to extend descriptions with other
    attributes unique to a content provider. These extensions are not expected to be
    natively understood under VRSATILE, but may be used for pre-negotiated exchange of
    message attributes when needed.
    """
    model_config = ConfigDict(
        extra='allow',
    )
    type: Literal['Extension'] = Field('Extension', description='MUST be "Extension".')
    name: str = Field(..., description='A name for the Extension')
    value: Optional[Union[float, str, bool, Dict[str, Any], List[Any]]] = Field(
        None, description='Any primitive or structured object'
    )


class _Entity(BaseModel):
    """Entity is the root class of `core` classes model - those that have identifiers
    and other general metadata like labels, xrefs, urls, descriptions, etc. All core
    classes descend from and inherit its attributes.
    """

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE)."
    )
    label: Optional[str] = Field(
        None,
        description='A primary label for the entity.'
    )
    description: Optional[str] = Field(
        None,
        description='A free-text description of the entity.'
    )
    extensions: Optional[List[Extension]] = None


class Coding(BaseModel):
    """a concept codified by a terminology system."""

    label: Optional[str] = Field(
        None,
        description='A primary label for the coding.'
    )
    system: str = Field(
        ...,
        description='Identity of the terminology system.'
    )
    version: Optional[str] = Field(
        None,
        description='Version of the terminology system.'
    )
    code: Code = Field(
        ...,
        description='Symbol in syntax defined by the terminology system.'
    )

class Mapping(_Entity):
    """A mapping to a concept in a terminology system."""
    model_config = ConfigDict(
        use_enum_values=True
    )

    coding: Coding
    relation: Relation = Field(
        ...,
        description='A mapping relation between concepts as defined by the Simple Knowledge Organization System (SKOS).'
    )


class _MappableEntity(_Entity):
    """an Entity that is mappable to codings in other terminology systems."""

    mappings: Optional[List[Mapping]] = None


class _DomainEntity(_MappableEntity):
    """An Entity that is specific to a particular biomedical domain such as disease,
    therapeutics, or genes.
    """

    type: str
    aliases: Optional[List[str]] = Field(
        None,
        description='Aliases are alternate labels for a Domain Entity.'
    )


class Phenotype(_DomainEntity):
    """An observable characteristic or trait of an organism."""

    type: Literal['Phenotype'] = Field(
        'Phenotype',
        description='MUST be "Phenotype".'
    )


class Disease(_DomainEntity):
    """A particular abnormal condition that negatively affects the structure or function
    of all or part of an organism and is not immediately due to any external injury.
    """

    type: Literal['Disease'] = Field(
        'Disease',
        description='MUST be "Disease".'
    )


class TraitSet(_DomainEntity):
    """A set of phenotype and/or disease concepts that together constitute a condition."""

    type: Literal['TraitSet'] = Field(
        'TraitSet',
        description='MUST be "TraitSet".'
    )
    traits: List[Union[Disease, Phenotype]] = Field(
        ...,
        min_length=2
    )


class TherapeuticAction(_DomainEntity):
    """A therapeutic action taken that is intended to alter or stop a pathologic process."""

    type: Literal['TherapeuticAction'] = Field(
        'TherapeuticAction',
        description='MUST be "TherapeuticAction".'
    )


class TherapeuticAgent(_DomainEntity):
    """An administered therapeutic agent that is intended to alter or stop a pathologic process."""

    type: Literal['TherapeuticAgent'] = Field(
        'TherapeuticAgent',
        description='MUST be "TherapeuticAgent".'
    )


class TherapeuticSubstituteGroup(_DomainEntity):
    """A group of therapeutic procedures that may be treated as substitutes for one another."""

    type: Literal['TherapeuticSubstituteGroup'] = Field(
        'TherapeuticSubstituteGroup',
        description='MUST be "TherapeuticSubstituteGroup".'
    )
    substitutes: List[Union[TherapeuticAction, TherapeuticAgent]] = Field(
        ...,
        description='The individual therapeutic procedures that may be treated as substitutes.',
        min_length=2
    )


class CombinationTherapy(_DomainEntity):
    """A therapeutic procedure that involves multiple different therapeutic procedures
    performed in combination.
    """

    type: Literal['CombinationTherapy'] = Field(
        'CombinationTherapy',
        description='MUST be "CombinationTherapy".'
    )
    components: List[Union[
        TherapeuticSubstituteGroup,
        TherapeuticAction,
        TherapeuticAgent
    ]] = Field(
        ...,
        description='The individual therapeutic procedure components that constitute the combination therapy.',
        min_length=2
    )


class Gene(_DomainEntity):
    """A basic physical and functional unit of heredity."""

    type: Literal['Gene'] = Field(
        'Gene',
        description='MUST be "Gene".'
    )


class _ValueObject(_Entity):
    """A contextual value whose equality is based on value, not identity. All VRS Value
    Objects may have computed digests from the VRS Computed Identifier algorithm.
    See https://en.wikipedia.org/wiki/Value_object for more on Value Objects.
    """

    digest: Optional[constr(pattern=r'[0-9A-Za-z_\-]{32}')] = Field(
        None,
        description='A sha512t24u digest created using the VRS Computed Identifier algorithm.',
    )


class _Ga4ghIdentifiableObject(_ValueObject):
    """A contextual value object for which a GA4GH computed identifier can be created."""

    type: str

    class ga4gh:
        identifiable = True
        prefix: str
        keys: List[str]


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
        description='An inclusive range of values bounded by one or more integers.',
        max_length=2,
        min_length=2,
    )


class Residue(RootModel):
    root: constr(pattern=r'[A-Z*\-]') = Field(
        ...,
        description='A character representing a specific residue (i.e., molecular species) or groupings of these ("ambiguity codes"), using [one-letter IUPAC abbreviations](https://en.wikipedia.org/wiki/International_Union_of_Pure_and_Applied_Chemistry#Amino_acid_and_nucleotide_base_codes) for nucleic acids and amino acids.',
    )


class SequenceString(RootModel):
    root: constr(pattern=r'^[A-Z*\-]*$') = Field(
        ...,
        description='A character string of Residues that represents a biological sequence using the conventional sequence order (5’-to-3’ for nucleic acid sequences, and amino-to-carboxyl for amino acid sequences). IUPAC ambiguity codes are permitted in Sequence Strings.',
    )


class SequenceReference(_ValueObject):
    model_config = ConfigDict(
        use_enum_values=True
    )

    type: Literal['SequenceReference'] = Field('SequenceReference', description='MUST be "SequenceReference"')
    refgetAccession: constr(pattern=r'SQ.[0-9A-Za-z_\-]{32}') = Field(
        ...,
        description='A `GA4GH RefGet <http://samtools.github.io/hts-specs/refget.html>` identifier for the referenced sequence, using the sha512t24u digest.',
    )
    residueAlphabet: Optional[ResidueAlphabet] = None

    class ga4gh:
        assigned: bool = Field(
            True,
            description='This special property indicates that the `digest` field follows an alternate convention and is expected to have the value assigned following that convention. For SequenceReference, it is expected the digest will be the refget accession value without the `SQ.` prefix.'
        )


class ReferenceLengthExpression(_ValueObject):
    """An expression of a length of a sequence from a repeating reference."""

    type: Literal['ReferenceLengthExpression'] = Field(
        'ReferenceLengthExpression', description='MUST be "ReferenceLengthExpression"'
    )
    length: Union[Range, int] = Field(
        ..., description='The number of residues in the expressed sequence.'
    )
    sequence: Optional[SequenceString] = Field(
        None, description='the Sequence encoded by the Reference Length Expression.'
    )
    repeatSubunitLength: int = Field(
        None, description='The number of residues in the repeat subunit.'
    )


class LiteralSequenceExpression(_ValueObject):
    """An explicit expression of a Sequence."""

    type: Literal['LiteralSequenceExpression'] = Field(
        'LiteralSequenceExpression', description='MUST be "LiteralSequenceExpression"'
    )
    sequence: SequenceString = Field(..., description='the literal sequence')


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

    class ga4gh:
        identifiable = True
        prefix = 'SL'
        keys = [
            'type',
            'start',
            'end',
            'sequenceReference'
        ]


class Allele(_Ga4ghIdentifiableObject):

    type: Literal['Allele'] = Field('Allele', description='MUST be "Allele"')
    location: Union[IRI, SequenceLocation] = Field(
        ..., description='The location of the Allele'
    )
    state: Union[LiteralSequenceExpression, ReferenceLengthExpression] = Field(
        ..., description='An expression of the sequence state'
    )

    class ga4gh:
        identifiable = True
        prefix = 'VA'
        keys = [
            'location',
            'state',
            'type'
        ]


class Haplotype(_Ga4ghIdentifiableObject):
    """A set of non-overlapping Allele members that co-occur on the same molecule."""

    type: Literal['Haplotype'] = Field('Haplotype', description='MUST be "Haplotype"')
    # TODO members temporarily typed as List instead of Set
    members: List[Union[Allele, IRI]] = Field(
        ...,
        description='A list of Alleles (or IRI references to `Alleles`) that comprise a Haplotype. Since each `Haplotype` member MUST be an `Allele`, and all members MUST share a common `SequenceReference`, implementations MAY use a compact representation of Haplotype that omits type and `SequenceReference` information in individual Haplotype members. Implementations MUST transform compact `Allele` representations into an `Allele` when computing GA4GH identifiers.',
        min_length=2,
    )

    class ga4gh:
        identifiable = True
        prefix = 'HT'
        keys = [
            'members',
            'type'
        ]


class _CopyNumber(_Ga4ghIdentifiableObject):
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

    class ga4gh:
        identifiable = True
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

    class ga4gh:
        identifiable = True
        prefix = 'CX'
        keys = [
            'copyChange',
            'location',
            'type'
        ]


class GenotypeMember(_ValueObject):
    """A class for expressing the count of a specific `MolecularVariation` present
    in-trans at a genomic locus represented by a `Genotype`.
    """

    type: Literal['GenotypeMember'] = Field('GenotypeMember', description='MUST be "GenotypeMember".')
    count: Union[Range, int] = Field(
        ..., description='The number of copies of the `variation` at a Genotype locus.'
    )
    variation: Union[Allele, Haplotype] = Field(
        ..., description='A MolecularVariation at a Genotype locus.'
    )

    class ga4gh:
        keys = [
            'type',
            'count',
            'variation'
        ]


class MolecularVariation(RootModel):
    """A variation on a contiguous molecule."""

    root: Union[Allele, Haplotype] = Field(
        ..., description='A variation on a contiguous molecule.', discriminator='type'
    )


class Genotype(_Ga4ghIdentifiableObject):
    """A quantified set of _in-trans_ `MolecularVariation` at a genomic locus."""

    type: Literal['Genotype'] = Field(
        'Genotype',
        description='MUST be "Genotype"'
    )
    # TODO members temporarily typed as List instead of Set + validate unique items
    members: List[GenotypeMember] = Field(
        ...,
        description='Each GenotypeMember in `members` describes a MolecularVariation and the count of that variation at the locus.',
        min_length=1,
    )
    count: Union[Range, int] = Field(
        ...,
        description='The total number of copies of all MolecularVariation at this locus, MUST be greater than or equal to the sum of GenotypeMember copy counts. If greater than the total counts, this implies additional MolecularVariation that are expected to exist but are not explicitly indicated.',
    )

    class ga4gh:
        identifiable = True
        prefix = 'GT'
        keys = [
            'count',
            'members',
            'type'
        ]


class SequenceExpression(RootModel):
    root: Union[LiteralSequenceExpression, ReferenceLengthExpression] = Field(
        ..., description='An expression describing a Sequence.', discriminator='type'
    )


class Location(RootModel):
    root: SequenceLocation = Field(
        ...,
        description='A contiguous segment of a biological sequence.',
        discriminator='type',
    )


class Condition(RootModel):
    """A disease or other medical disorder."""

    root: Union[Disease, Phenotype, TraitSet] = Field(
        ...,
        description='A disease or other medical disorder.',
        discriminator='type',
    )


class TherapeuticProcedure(RootModel):
    """An action or administration of therapeutic agents to produce an effect  that is
    intended to alter or stop a pathologic process.
    """

    root: Union[CombinationTherapy, TherapeuticAction, TherapeuticAgent, TherapeuticSubstituteGroup] = Field(
        ...,
        description='An action or administration of therapeutic agents to produce an effect that is intended to alter or stop a pathologic process.',
        discriminator='type',
    )


class Variation(RootModel):
    root: Union[Allele, CopyNumberChange, CopyNumberCount, Genotype, Haplotype] = Field(
        ...,
        description='A representation of the state of one or more biomolecules.',
        discriminator='type',
    )


class SystemicVariation(RootModel):
    """A Variation of multiple molecules in the context of a system, e.g. a genome,
    sample, or homologous chromosomes.
    """

    root: Union[CopyNumberChange, CopyNumberCount, Genotype] = Field(
        ...,
        description='A Variation of multiple molecules in the context of a system, e.g. a genome, sample, or homologous chromosomes.',
        discriminator='type',
    )


# At end so classes exist
(
    reffable_classes,
    union_reffable_classes,
    class_refatt_map,
    class_keys
) = pydantic_class_refatt_map()
