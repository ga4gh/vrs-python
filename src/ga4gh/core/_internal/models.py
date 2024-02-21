"""
**This module should not be imported directly.**

Instead, users should use one of the following:

  * `from ga4gh.core import core_models`, and refer to models with the
    abbreviated name, e.g., `core_models.Gene` (recommended)

  * `import ga4gh.core`, and refer to models using the fully-qualified
    module name, e.g., `ga4gh.core.core_models.Gene`
"""
from typing import Any, Dict, List, Literal, Optional, Union
from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, RootModel, constr, model_serializer
from ga4gh.core import GA4GH_IR_REGEXP


class Relation(Enum):
    """A mapping relation between concepts as defined by the Simple Knowledge
    Organization System (SKOS).
    """

    CLOSE_MATCH = 'closeMatch'
    EXACT_MATCH = 'exactMatch'
    BROAD_MATCH = 'broadMatch'
    NARROW_MATCH = 'narrowMatch'
    RELATED_MATCH = 'relatedMatch'


class Code(RootModel):
    """Indicates that the value is taken from a set of controlled strings defined
    elsewhere. Technically, a code is restricted to a string which has at least one
    character and no leading or  trailing whitespace, and where there is no whitespace
    other than single spaces in the contents."""

    root: constr(pattern=r'\S+( \S+)*') = Field(
        ...,
        json_schema_extra={
            'description': 'Indicates that the value is taken from a set of controlled strings defined elsewhere. Technically, a code is restricted to a string which has at least one character and no leading or  trailing whitespace, and where there is no whitespace other than single spaces in the contents.',
            'example': 'ENSG00000139618',
        }
    )


class IRI(RootModel):

    def __hash__(self):
        return self.root.__hash__()

    @model_serializer(when_used='json')
    def ga4gh_serialize(self):
        m = GA4GH_IR_REGEXP.match(self.root)
        if m is not None:
            return m['digest']
        return self.root

    root: str = Field(
        ...,
        json_schema_extra={'description': 'An IRI Reference (either an IRI or a relative-reference), according to `RFC3986 section 4.1  <https://datatracker.ietf.org/doc/html/rfc3986#section-4.1>` and `RFC3987 section 2.1 <https://datatracker.ietf.org/doc/html/rfc3987#section-2.1>`. MAY be a JSON Pointer as an IRI fragment, as  described by `RFC6901 section 6 <https://datatracker.ietf.org/doc/html/rfc6901#section-6>`.',
        }
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


class Condition(RootModel):
    """A disease or other medical disorder."""

    root: Union[Disease, Phenotype, TraitSet] = Field(
        ...,
        json_schema_extra={'description': 'A disease or other medical disorder.'},
        discriminator='type',
    )


class TherapeuticProcedure(RootModel):
    """An action or administration of therapeutic agents to produce an effect  that is
    intended to alter or stop a pathologic process.
    """

    root: Union[CombinationTherapy, TherapeuticAction, TherapeuticAgent, TherapeuticSubstituteGroup] = Field(
        ...,
        json_schema_extra={'description': 'An action or administration of therapeutic agents to produce an effect that is intended to alter or stop a pathologic process.'},
        discriminator='type',
    )


class Gene(_DomainEntity):
    """A basic physical and functional unit of heredity."""

    type: Literal['Gene'] = Field(
        'Gene',
        description='MUST be "Gene".'
    )
