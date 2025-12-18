"""GKS Core Class Definitions"""

from __future__ import annotations

from abc import ABC
from enum import StrEnum
from typing import Annotated, Any, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    RootModel,
    StringConstraints,
    model_validator,
)
from typing_extensions import Self

from ga4gh.core.identifiers import GA4GH_IR_REGEXP


class BaseModelForbidExtra(BaseModel):
    """Base Pydantic model class with extra attributes forbidden."""

    model_config = ConfigDict(extra="forbid")


class Relation(StrEnum):
    """A mapping relation between concepts as defined by the Simple Knowledge
    Organization System (SKOS).
    """

    CLOSE_MATCH = "closeMatch"
    EXACT_MATCH = "exactMatch"
    BROAD_MATCH = "broadMatch"
    NARROW_MATCH = "narrowMatch"
    RELATED_MATCH = "relatedMatch"


class MembershipOperator(StrEnum):
    """The logical relationship between concepts in the set, in the context of some
    knowledge reported about them. The value 'AND' indicates that the concepts are
    dependent and occur together in this context - i.e. the reported assertion is not
    necessarily true for each concept on its own - only in combination with the
    other(s). The value 'OR' indicates that each concept applies independently in this
    context - i.e. the reported assertion is necessarily true for each concept on its
    own, independent of the presence of the other(s).
    """

    AND = "AND"
    OR = "OR"


#########################################
# Primitive data types
#########################################


class code(RootModel):  # noqa: N801
    """Indicates that the value is taken from a set of controlled strings defined
    elsewhere. Technically, a code is restricted to a string which has at least one
    character and no leading or trailing whitespace, and where there is no whitespace
    other than single spaces in the contents.
    """

    root: Annotated[str, StringConstraints(pattern=r"\S+( \S+)*")] = Field(
        ...,
        json_schema_extra={
            "description": "Indicates that the value is taken from a set of controlled strings defined elsewhere. Technically, a code is restricted to a string which has at least one character and no leading or  trailing whitespace, and where there is no whitespace other than single spaces in the contents.",
            "example": "ENSG00000139618",
        },
    )


class iriReference(RootModel):  # noqa: N801
    """An IRI Reference (either an IRI or a relative-reference), according to `RFC3986
    section 4.1 <https://datatracker.ietf.org/doc/html/rfc3986#section-4.1>`_ and
    `RFC3987 section 2.1 <https://datatracker.ietf.org/doc/html/rfc3987#section-2.1>`_.
    MAY be a JSON Pointer as an IRI fragment, as described by `RFC6901 section 6
    <https://datatracker.ietf.org/doc/html/rfc6901#section-6>`_.
    """

    def __hash__(self) -> int:  # noqa: D105
        return self.root.__hash__()

    def ga4gh_serialize(self) -> str:  # noqa: D102
        m = GA4GH_IR_REGEXP.match(self.root)
        if m is not None:
            return m["digest"]
        return self.root

    root: str = Field(
        ...,
        json_schema_extra={
            "description": "An IRI Reference (either an IRI or a relative-reference), according to `RFC3986 section 4.1 <https://datatracker.ietf.org/doc/html/rfc3986#section-4.1>`_ and `RFC3987 section 2.1 <https://datatracker.ietf.org/doc/html/rfc3987#section-2.1>`_. MAY be a JSON Pointer as an IRI fragment, as described by `RFC6901 section 6 <https://datatracker.ietf.org/doc/html/rfc6901#section-6>`_.",
        },
    )


#########################################
# Abstract core classes
#########################################


class Entity(BaseModel, ABC):
    """Anything that exists, has existed, or will exist.

    Abstract base class to be extended by other classes. Do NOT instantiate directly.
    """

    id: str | None = Field(
        default=None,
        description="The 'logical' identifier of the Entity in the system of record, e.g. a UUID.  This 'id' is unique within a given system, but may or may not be globally unique outside the system. It is used within a system to reference an object from another.",
    )
    type: str = Field(
        ...,
        description="The name of the class that is instantiated by a data object representing the Entity.",
    )
    name: str | None = Field(default=None, description="A primary name for the entity.")
    description: str | None = Field(
        default=None, description="A free-text description of the Entity."
    )
    aliases: list[str] | None = Field(
        default=None, description="Alternative name(s) for the Entity."
    )
    extensions: list[Extension] | None = Field(
        default=None,
        description="A list of extensions to the Entity, that allow for capture of information not directly supported by elements defined in the model.",
    )


class Element(BaseModel, ABC):
    """The base definition for all identifiable data objects.

    Abstract base class to be extended by other classes. Do NOT instantiate directly.
    """

    id: str | None = Field(
        default=None,
        description="The 'logical' identifier of the data element in the system of record, e.g. a UUID.  This 'id' is unique within a given system, but may or may not be globally unique outside the system. It is used within a system to reference an object from another.",
    )
    extensions: list[Extension] | None = Field(
        default=None,
        description="A list of extensions to the Entity, that allow for capture of information not directly supported by elements defined in the model.",
    )


#########################################
# General-purpose data classes
#########################################


class Coding(Element, BaseModelForbidExtra):
    """A structured representation of a code for a defined concept in a terminology or
    code system.
    """

    name: str | None = Field(
        default=None,
        description="The human-readable name for the coded concept, as defined by the code system.",
    )
    system: str = Field(
        ...,
        description="The terminology/code system that defined the code. May be reported as a free-text name (e.g. 'Sequence Ontology'), but it is preferable to provide a uri/url for the system.",
    )
    systemVersion: str | None = Field(  # noqa: N815
        default=None,
        description="Version of the terminology or code system that provided the code.",
    )
    code: code  # Cannot use Field due to PydanticUserError: field name and type annotation must not clash.
    iris: list[iriReference] | None = Field(
        default=None,
        description="A list of IRIs that are associated with the coding. This can be used to provide additional context or to link to additional information about the concept.",
    )


class ConceptMapping(Element, BaseModelForbidExtra):
    """A mapping to a concept in a terminology or code system."""

    model_config = ConfigDict(use_enum_values=True)

    coding: Coding = Field(
        ...,
        description="A structured representation of a code for a defined concept in a terminology or code system.",
    )
    relation: Relation = Field(
        ...,
        description="A mapping relation between concepts as defined by the Simple Knowledge Organization System (SKOS).",
    )


class ConceptSet(Element, BaseModelForbidExtra):
    """A set of concepts that may be considered as dependent (occurring together), or
    independent (existing separately) in the context of some knowledge reported about
    them, as indicated by a set membership operator. e.g. a set of independent molecular
    consequences that both result from the presence of a particular genetic variant
    (membership operator = OR).
    """

    type: Literal["ConceptSet"] = Field(
        default="ConceptSet",
        description='MUST be "ConceptSet"',
    )
    concepts: list[MappableConcept] | list[ConceptSet] = Field(
        ...,
        description="A list of concepts that are dependent (occurring together), or independent (existing separately), depending on the",
        min_length=2,
    )
    membershipOperator: MembershipOperator = Field(  # noqa: N815
        ...,
        description="The logical relationship between concepts in the set, in the context of some knowledge reported about them. The value 'AND' indicates that the concepts are dependent and occur together in this context - i.e. the reported assertion is not necessarily true for each concept on its own - only in combination with the other(s). The value 'OR' indicates that each concept applies independently in this context - i.e. the reported assertion is necessarily true for each concept on its own, independent of the presence of the other(s).",
    )


class Extension(Element, BaseModelForbidExtra):
    """The Extension class provides entities with a means to include additional
    attributes that are outside of the specified standard but needed by a given content
    provider or system implementer. These extensions are not expected to be natively
    understood, but may be used for pre-negotiated exchange of message attributes
    between systems.
    """

    name: str = Field(
        ...,
        description="A name for the Extension. Should be indicative of its meaning and/or the type of information it value represents.",
    )
    value: float | str | bool | dict[str, Any] | list[Any] | None = Field(
        ...,
        description="The value of the Extension - can be any primitive or structured object",
    )
    description: str | None = Field(
        default=None,
        description="A description of the meaning or utility of the Extension, to explain the type of information it is meant to hold.",
    )


class MappableConcept(Element, BaseModelForbidExtra):
    """A concept based on a primaryCoding and/or name that may be mapped to one or more other `Codings`."""

    conceptType: str | None = Field(  # noqa: N815
        default=None,
        description="A term indicating the type of concept being represented by the MappableConcept.",
    )
    name: str | None = Field(
        default=None, description="A primary name for the concept."
    )
    primaryCoding: Coding | None = Field(  # noqa: N815
        default=None,
        description="A primary coding for the concept.",
    )
    mappings: list[ConceptMapping] | None = Field(
        default=None,
        description="A list of mappings to concepts in terminologies or code systems. Each mapping should include a coding and a relation.",
    )

    @model_validator(mode="after")
    def require_name_or_primary_coding(self) -> Self:
        """Ensure that ``name`` or ``primaryCoding`` is provided"""
        if self.primaryCoding is None and self.name is None:
            err_msg = "One of `name` or `primaryCoding` must be provided."
            raise ValueError(err_msg)
        return self


Element.model_rebuild()
Entity.model_rebuild()
