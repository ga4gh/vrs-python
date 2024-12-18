"""GKS Core Class Definitions"""
from __future__ import annotations

from abc import ABC
from datetime import date as datetime_date, datetime as datetime_datetime
from typing import Any, Dict, Annotated, Optional, Union, List
from enum import Enum

from pydantic import BaseModel, Field, RootModel, StringConstraints, ConfigDict, model_validator

from ga4gh.core import GA4GH_IR_REGEXP


class Relation(str, Enum):
    """A mapping relation between concepts as defined by the Simple Knowledge
    Organization System (SKOS).
    """

    CLOSE_MATCH = 'closeMatch'
    EXACT_MATCH = 'exactMatch'
    BROAD_MATCH = 'broadMatch'
    NARROW_MATCH = 'narrowMatch'
    RELATED_MATCH = 'relatedMatch'


#########################################
# Primitive data types
#########################################


class code(RootModel):
    """Indicates that the value is taken from a set of controlled strings defined
    elsewhere. Technically, a code is restricted to a string which has at least one
    character and no leading or trailing whitespace, and where there is no whitespace
    other than single spaces in the contents."""

    root: Annotated[str, StringConstraints(pattern=r'\S+( \S+)*')] = Field(
        ...,
        json_schema_extra={
            'description': 'Indicates that the value is taken from a set of controlled strings defined elsewhere. Technically, a code is restricted to a string which has at least one character and no leading or  trailing whitespace, and where there is no whitespace other than single spaces in the contents.',
            'example': 'ENSG00000139618',
        }
    )


class iriReference(RootModel):
    """An IRI Reference (either an IRI or a relative-reference), according to `RFC3986
    section 4.1 <https://datatracker.ietf.org/doc/html/rfc3986#section-4.1>`_ and
    `RFC3987 section 2.1 <https://datatracker.ietf.org/doc/html/rfc3987#section-2.1>`_.
    MAY be a JSON Pointer as an IRI fragment, as described by `RFC6901 section 6
    <https://datatracker.ietf.org/doc/html/rfc6901#section-6>`_.
    """

    def __hash__(self):
        return self.root.__hash__()

    def ga4gh_serialize(self):
        m = GA4GH_IR_REGEXP.match(self.root)
        if m is not None:
            return m['digest']
        return self.root

    root: str = Field(
        ...,
        json_schema_extra={'description': 'An IRI Reference (either an IRI or a relative-reference), according to `RFC3986 section 4.1 <https://datatracker.ietf.org/doc/html/rfc3986#section-4.1>`_ and `RFC3987 section 2.1 <https://datatracker.ietf.org/doc/html/rfc3987#section-2.1>`_. MAY be a JSON Pointer as an IRI fragment, as described by `RFC6901 section 6 <https://datatracker.ietf.org/doc/html/rfc6901#section-6>`_.',
        }
    )


class date(RootModel):
    """A string is valid against this format if it represents a date in the following format: YYYY-MM-DD."""

    root: datetime_date = Field(..., json_schema_extra={"description": "A string is valid against this format if it represents a date in the following format: YYYY-MM-DD."})


class datetime(RootModel):
    """A string is valid against this format if it represents a date-time in the
    following format: YYYY:MM::DDThh:mm:ss.sTZD.."""

    root: datetime_datetime = Field(..., json_schema_extra={"description": "A string is valid against this format if it represents a date-time in the following format: YYYY:MM::DDThh:mm:ss.sTZD."})


#########################################
# Abstract core classes
#########################################


class Entity(BaseModel, ABC):
    """Anything that exists, has existed, or will exist.

    Abstract base class to be extended by other classes. Do NOT instantiate directly.
    """

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the Entity in the system of record, e.g. a UUID.  This 'id' is unique within a given system, but may or may not be globally unique outside the system. It is used within a system to reference an object from another."
    )
    type: str = Field(..., description="The name of the class that is instantiated by a data object representing the Entity.")
    label: Optional[str] = Field(
        None,
        description='A primary name for the entity.'
    )
    description: Optional[str] = Field(
        None,
        description='A free-text description of the Entity.'
    )
    alternativeLabels: Optional[List[str]] = Field(None, description="Alternative name(s) for the Entity.")
    extensions: Optional[List[Extension]] = Field(None, description="A list of extensions to the Entity, that allow for capture of information not directly supported by elements defined in the model.")


class Element(BaseModel, ABC):
    """The base definition for all identifiable data objects.

    Abstract base class to be extended by other classes. Do NOT instantiate directly.
    """

    id: Optional[str] = Field(None, description="The 'logical' identifier of the data element in the system of record, e.g. a UUID.  This 'id' is unique within a given system, but may or may not be globally unique outside the system. It is used within a system to reference an object from another.")
    extensions: Optional[List[Extension]] = Field(None, description="A list of extensions to the Entity, that allow for capture of information not directly supported by elements defined in the model.")


#########################################
# General-purpose data classes
#########################################

class Coding(Element):
    """A structured representation of a code for a defined concept in a terminology or
    code system.
    """

    label: Optional[str] = Field(
        None,
        description='The human-readable name for the coded concept, as defined by the code system.'
    )
    system: str = Field(
        ...,
        description="The terminology/code system that defined the code. May be reported as a free-text name (e.g. 'Sequence Ontology'), but it is preferable to provide a uri/url for the system. When the 'code' is reported as a CURIE, the 'system' should be reported as the uri that the CURIE's prefix expands to (e.g. 'http://purl.obofoundry.org/so.owl/' for the Sequence Ontology)."
    )
    systemVersion: Optional[str] = Field(
        None,
        description='Version of the terminology or code system that provided the code.'
    )
    code: "code"  # Cannot use Field due to PydanticUserError: field name and type annotation must not clash.


class ConceptMapping(Element):
    """A mapping to a concept in a terminology or code system."""

    model_config = ConfigDict(use_enum_values=True)

    coding: Coding = Field(..., description="A structured representation of a code for a defined concept in a terminology or code system.")
    relation: Relation = Field(..., description="A mapping relation between concepts as defined by the Simple Knowledge Organization System (SKOS).")


class Extension(Element):
    """The Extension class provides entities with a means to include additional
    attributes that are outside of the specified standard but needed by a given content
    provider or system implementer. These extensions are not expected to be natively
    understood, but may be used for pre-negotiated exchange of message attributes
    between systems.
    """

    name: str = Field(..., description='A name for the Extension. Should be indicative of its meaning and/or the type of information it value represents.')
    value: Optional[Union[float, str, bool, Dict[str, Any], List[Any]]] = Field(
        ..., description='The value of the Extension - can be any primitive or structured object'
    )
    description: Optional[str] = Field(None, description="A description of the meaning or utility of the Extension, to explain the type of information it is meant to hold.")


class MappableConcept(Element):
    """A concept label that may be mapped to one or more :ref:`Codings <Coding>`."""

    conceptType: Optional[str] = Field(None, description="A term indicating the type of concept being represented by the MappableConcept.")
    label: Optional[str] = Field(None, description="A primary name for the concept.")
    primaryCode: Optional[code] = Field(None, description="A primary code for the concept that is used to identify the concept in a terminology or code system. If there is a public code system for the primaryCode then it should also be specified in the mappings array with a relation of 'exactMatch'. This attribute is provided to both allow a more technical code to be used when a public Coding with a system is not available as well as when it is available but should be identified as the primary code.")
    mappings: Optional[List[ConceptMapping]] = Field(None, description="A list of mappings to concepts in terminologies or code systems. Each mapping should include a coding and a relation.")

    @model_validator(mode="after")
    def require_label_or_primary_code(cls, v):
        """Ensure that ``label`` or ``primaryCode`` is provided"""
        if v.primaryCode is None and v.label is None:
            err_msg = "`label` or `primaryCode` must be provided."
            raise ValueError(err_msg)
        return v

    def ga4gh_serialize(self):
        return self.primaryCode.root

Element.model_rebuild()
Entity.model_rebuild()
Extension.model_rebuild()
