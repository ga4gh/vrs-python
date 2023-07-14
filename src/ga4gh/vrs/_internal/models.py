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
datamodel-codegen --input submodules/vrs/schema/merged.json --input-file-type jsonschema --output models_merged2.py

"""
from typing import Any, Dict, List, Optional, Union, Literal
from enum import Enum
import inspect
import logging
import sys
import os
import typing
import pkg_resources
from pydantic import BaseModel, Extra, Field, constr
_logger = logging.getLogger(__name__)


# specify VRSATILE_SCHEMA_DIR to use a schema other than the one embedded
# in vrs-python
schema_dir = os.environ.get("VRSATILE_SCHEMA_DIR", pkg_resources.resource_filename(__name__, "data/schemas/vrsatile"))
schema_path = schema_dir + "/merged.json"


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
    if ('__origin__' in t.__dict__):
        if t.__origin__ == typing.Literal:
            return list(t.__args__)
        elif (t.__origin__ == typing.Union
              or issubclass(t.__origin__,  typing.List)):
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
        lambda c: ('id' in c.__fields__
                   and hasattr(c, 'ga4gh_digest')),
        model_classes
    ))
    # Types reffable because they are a union of reffable types
    union_reffable_classes = list(filter(
        lambda c: ('__root__' in c.__fields__
                   and overlaps(reffable_classes,
                                flatten_type(c.__fields__['__root__'].type_))),
        model_classes
    ))
    reffable_fields = {}
    # Find any field whose type is a subclass of a reffable type,
    # or which is a typing.List that includes a reffable type.
    # Interestingly, ModelField.type_ is the member type for List types, not type List itself.
    for model_class in model_classes:
        # if model_class in union_reffable_classes:
        #     continue  # These don't have fields other than __root__
        fields = model_class.__fields__
        class_reffable_fields = []
        for fieldname, field in fields.items():
            if fieldname == '__root__':
                continue
            field_type = field.type_  # a typing or normal annotation like str
            # types can be raw class type annotation (int, str, dict, etc)
            # typing.Literal, typing.Union, typing.Optional
            # Use flatten_type to simplify these
            if any([rc in flatten_type(field_type)
                    for rc in (reffable_classes
                               + union_reffable_classes)]):
                class_reffable_fields.append(fieldname)
        if len(class_reffable_fields) > 0:
            reffable_fields[model_class.__name__] = class_reffable_fields
    return reffable_fields


class Extension(BaseModel):
    class Config:
        extra = Extra.forbid

    type: str = Field('Extension', const=True, description='MUST be "Extension".')
    name: str = Field(..., description='A name for the Extension')
    value: Optional[Union[float, str, bool, Dict[str, Any], List[Any]]] = Field(
        None, description='Any primitive or structured object'
    )


class MappingClass(Enum):
    closeMatch = 'closeMatch'
    exactMatch = 'exactMatch'
    broadMatch = 'broadMatch'
    narrowMatch = 'narrowMatch'
    relatedMatch = 'relatedMatch'


class Code(BaseModel):
    __root__: constr(regex=r'\S+( \S+)*') = Field(
        ...,
        description='Indicates that the value is taken from a set of controlled strings defined elsewhere. Technically, a code is restricted to a string which has at least one character and no leading or  trailing whitespace, and where there is no whitespace other than single spaces in the contents.',
        example='ENSG00000139618',
    )


class IRI(BaseModel):
    __root__: str = Field(
        ...,
        description='An IRI Reference (either an IRI or a relative-reference), according to `RFC3986 section 4.1  <https://datatracker.ietf.org/doc/html/rfc3986#section-4.1>` and `RFC3987 section 2.1 <https://datatracker.ietf.org/doc/html/rfc3987#section-2.1>`. MAY be a JSON Pointer as an IRI fragment, as  described by `RFC6901 section 6 <https://datatracker.ietf.org/doc/html/rfc6901#section-6>`.',
    )


class CopyChange(Enum):
    efo_0030069 = 'efo:0030069'
    efo_0020073 = 'efo:0020073'
    efo_0030068 = 'efo:0030068'
    efo_0030067 = 'efo:0030067'
    efo_0030064 = 'efo:0030064'
    efo_0030070 = 'efo:0030070'
    efo_0030071 = 'efo:0030071'
    efo_0030072 = 'efo:0030072'


class ResidueAlphabet(Enum):
    amino_acid = 'amino acid'
    nucleic_acid = 'nucleic acid'


class SequenceReference(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE).",
    )
    label: Optional[str] = None
    extensions: Optional[List[Extension]] = None
    type: str = Field(
        'SequenceReference', const=True, description='MUST be "SequenceReference"'
    )
    digest: Optional[constr(regex=r'[0-9A-Za-z_\-]{32}')] = Field(
        None,
        description='A sha512t24u digest created using the VRS Computed Identifier algorithm.',
    )
    refgetAccession: constr(regex=r'SQ.[0-9A-Za-z_\-]{32}') = Field(
        ...,
        description='A `GA4GH RefGet <http://samtools.github.io/hts-specs/refget.html>` identifier for the referenced sequence, using the sha512t24u digest.',
    )
    residueAlphabet: Optional[ResidueAlphabet] = None

    class ga4gh_digest:
        prefix = 'SQR'
        keys = [
            'refgetAccession',
            'type'
        ]


class Range(BaseModel):
    __root__: List[Optional[int]] = Field(
        ...,
        description='An inclusive range of values bounded by one or more integers.',
        max_items=2,
        min_items=2,
    )


class Residue(BaseModel):
    __root__: constr(regex=r'[A-Z*\-]') = Field(
        ...,
        description='A character representing a specific residue (i.e., molecular species) or groupings of these ("ambiguity codes"), using [one-letter IUPAC abbreviations](https://en.wikipedia.org/wiki/International_Union_of_Pure_and_Applied_Chemistry#Amino_acid_and_nucleotide_base_codes) for nucleic acids and amino acids.',
    )


class SequenceString(BaseModel):
    __root__: constr(regex=r'^[A-Z*\-]*$') = Field(
        ...,
        description='A character string of Residues that represents a biological sequence using the conventional sequence order (5’-to-3’ for nucleic acid sequences, and amino-to-carboxyl for amino acid sequences). IUPAC ambiguity codes are permitted in Sequence Strings.',
    )


class Mapping(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE).",
    )
    label: Optional[str] = None
    extensions: Optional[List[Extension]] = None
    system: str = Field(..., description='Identity of the terminology system.')
    version: Optional[str] = Field(
        None, description='Version of the terminology system.'
    )
    code: Code = Field(
        ..., description='Symbol in syntax defined by the terminology system.'
    )
    mapping: MappingClass = Field(
        ...,
        description='A mapping between concepts as defined by the Simple Knowledge Organization System (SKOS).',
    )


class SequenceLocation(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE).",
    )
    label: Optional[str] = None
    extensions: Optional[List[Extension]] = None
    type: str = Field(
        'SequenceLocation', const=True, description='MUST be "SequenceLocation"'
    )
    digest: Optional[constr(regex=r'[0-9A-Za-z_\-]{32}')] = Field(
        None,
        description='A sha512t24u digest created using the VRS Computed Identifier algorithm.',
    )
    sequence: Optional[Union[SequenceReference, IRI]] = Field(
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

    class ga4gh_digest:
        prefix = 'SL',
        keys = [
            'count',
            'members',
            'type'
        ]



class ReferenceLengthExpression(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE).",
    )
    label: Optional[str] = None
    extensions: Optional[List[Extension]] = None
    type: Literal['ReferenceLengthExpression'] = Field(
        'ReferenceLengthExpression',
        const=True,
        description='MUST be "ReferenceLengthExpression"',
    )
    length: Union[Range, int] = Field(
        ..., description='The number of residues in the expressed sequence.'
    )
    sequence: Optional[SequenceString] = Field(
        None, description='the Sequence encoded by the Reference Length Expression.'
    )
    repeatSubunitLength: Optional[int] = Field(
        None, description='The number of residues in the repeat subunit.'
    )


class LiteralSequenceExpression(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE).",
    )
    label: Optional[str] = None
    extensions: Optional[List[Extension]] = None
    type: Literal['LiteralSequenceExpression'] = Field(
        'LiteralSequenceExpression',
        const=True,
        description='MUST be "LiteralSequenceExpression"',
    )
    sequence: SequenceString = Field(..., description='the literal sequence')


class Allele(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE).",
    )
    label: Optional[str] = None
    extensions: Optional[List[Extension]] = None
    type: Literal['Allele'] = Field('Allele', const=True, description='MUST be "Allele"')
    digest: Optional[constr(regex=r'[0-9A-Za-z_\-]{32}')] = Field(
        None,
        description='A sha512t24u digest created using the VRS Computed Identifier algorithm.',
    )
    location: Union[SequenceLocation, IRI] = Field(
        ..., description='The location of the Allele'
    )
    state: Union[LiteralSequenceExpression, ReferenceLengthExpression] = Field(
        ..., description='An expression of the sequence state'
    )

    class ga4gh_digest:
        prefix = 'VA'
        keys = [
            'location',
            'state',
            'type'
        ]


class Haplotype(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE).",
    )
    label: Optional[str] = None
    extensions: Optional[List[Extension]] = None
    type: Literal['Haplotype'] = Field('Haplotype', const=True, description='MUST be "Haplotype"')
    digest: Optional[constr(regex=r'[0-9A-Za-z_\-]{32}')] = Field(
        None,
        description='A sha512t24u digest created using the VRS Computed Identifier algorithm.',
    )
    members: List[Union[Allele, IRI]] = Field(
        ...,
        description='List of Alleles, or references to Alleles, that comprise this Haplotype.',
        min_items=2,
        unique_items=True,
    )

    class ga4gh_digest:
        prefix = 'HT'
        keys = [
            'members',
            'type'
        ]



class CopyNumberCount(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE).",
    )
    label: Optional[str] = None
    extensions: Optional[List[Extension]] = None
    type: Literal['CopyNumberCount'] = Field(
        'CopyNumberCount', const=True, description='MUST be "CopyNumberCount"'
    )
    digest: Optional[constr(regex=r'[0-9A-Za-z_\-]{32}')] = Field(
        None,
        description='A sha512t24u digest created using the VRS Computed Identifier algorithm.',
    )
    subject: Union[SequenceLocation, IRI] = Field(
        ...,
        description='A location for which the number of systemic copies is described.',
    )
    copies: Union[Range, int] = Field(
        ..., description='The integral number of copies of the subject in a system'
    )

    class ga4gh_digest:
        prefix = 'CN'
        keys = [
            'copies',
            'subject',
            'type'
        ]


class CopyNumberChange(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE).",
    )
    label: Optional[str] = None
    extensions: Optional[List[Extension]] = None
    type: Literal['CopyNumberChange'] = Field(
        'CopyNumberChange', const=True, description='MUST be "CopyNumberChange"'
    )
    digest: Optional[constr(regex=r'[0-9A-Za-z_\-]{32}')] = Field(
        None,
        description='A sha512t24u digest created using the VRS Computed Identifier algorithm.',
    )
    subject: Union[SequenceLocation, IRI] = Field(
        ...,
        description='A location for which the number of systemic copies is described.',
    )
    copy_change: CopyChange = Field(
        ...,
        description='MUST be one of "efo:0030069" (complete genomic loss), "efo:0020073" (high-level loss), "efo:0030068" (low-level loss), "efo:0030067" (loss), "efo:0030064" (regional base ploidy), "efo:0030070" (gain), "efo:0030071" (low-level gain), "efo:0030072" (high-level gain).',
    )

    class ga4gh_digest:
        prefix = 'CX',
        keys = [
            'copy_change',
            'subject',
            'type'
        ]



class Location(BaseModel):
    __root__: SequenceLocation = Field(
        ...,
        description='A contiguous segment of a biological sequence.'
    )


class SequenceExpression(BaseModel):
    __root__: Union[LiteralSequenceExpression, ReferenceLengthExpression] = Field(
        ..., description='An expression describing a Sequence.', discriminator='type'
    )


class GenotypeMember(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE).",
    )
    label: Optional[str] = None
    extensions: Optional[List[Extension]] = None
    type: Literal['GenotypeMember'] = Field(
        'GenotypeMember', const=True, description='MUST be "GenotypeMember".'
    )
    count: Union[Range, int] = Field(
        ..., description='The number of copies of the `variation` at a Genotype locus.'
    )
    variation: Union[Allele, Haplotype] = Field(
        ..., description='A MolecularVariation at a Genotype locus.'
    )


class MolecularVariation(BaseModel):
    __root__: Union[Allele, Haplotype] = Field(
        ..., description='A variation on a contiguous molecule.', discriminator='type'
    )


class Genotype(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is  unique within a given system. The identified entity may have a different 'id' in a different  system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE).",
    )
    label: Optional[str] = None
    extensions: Optional[List[Extension]] = None
    type: Literal['Genotype'] = Field('Genotype', const=True, description='MUST be "Genotype"')
    digest: Optional[constr(regex=r'[0-9A-Za-z_\-]{32}')] = Field(
        None,
        description='A sha512t24u digest created using the VRS Computed Identifier algorithm.',
    )
    members: List[GenotypeMember] = Field(
        ...,
        description='Each GenotypeMember in `members` describes a MolecularVariation and the count of that variation at the locus.',
        min_items=1,
        unique_items=True,
    )
    count: Union[Range, int] = Field(
        ...,
        description='The total number of copies of all MolecularVariation at this locus, MUST be greater than or equal to the sum of GenotypeMember copy counts. If greater than the total counts, this implies additional MolecularVariation that are expected to exist but are not explicitly indicated.',
    )

    class ga4gh_digest:
        prefix = 'GT',
        keys = [
            'count',
            'members',
            'type'
        ]



class Variation(BaseModel):
    __root__: Union[
        Allele, CopyNumberChange, CopyNumberCount, Genotype, Haplotype
    ] = Field(
        ...,
        description='A representation of the state of one or more biomolecules.',
        discriminator='type',
    )


class SystemicVariation(BaseModel):
    __root__: Union[CopyNumberChange, CopyNumberCount, Genotype] = Field(
        ...,
        description='A Variation of multiple molecules in the context of a system, e.g. a genome, sample, or homologous chromosomes.',
        discriminator='type',
    )


# At end so classes exist
class_refatt_map = pydantic_class_refatt_map()
