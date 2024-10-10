"""GKS Common Library Data Type and Entity models"""
from __future__ import annotations

import datetime
import logging
from typing import Any, Dict, Annotated, Literal, Optional, Union, List
from enum import Enum

from pydantic import BaseModel, Field, RootModel, StringConstraints, ConfigDict, field_validator

from ga4gh.core import GA4GH_IR_REGEXP



class CoreImType(str, Enum):
    """Define Core Information Model Types"""

    AGENT = "Agent"
    CONTRIBUTION = "Contribution"
    DOCUMENT = "Document"
    METHOD = "Method"
    DATA_SET = "DataSet"
    EVIDENCE_LINE = "EvidenceLine"
    INFORMATION_ENTITY = "InformationEntity"
    STUDY_GROUP = "StudyGroup"



class Relation(str, Enum):
    """A mapping relation between concepts as defined by the Simple Knowledge
    Organization System (SKOS).
    """

    CLOSE_MATCH = 'closeMatch'
    EXACT_MATCH = 'exactMatch'
    BROAD_MATCH = 'broadMatch'
    NARROW_MATCH = 'narrowMatch'
    RELATED_MATCH = 'relatedMatch'


class AgentSubtype(str, Enum):
    """A specific type of agent the Agent object represents."""

    PERSON = "person"
    ORGANIZATION = "organization"
    SOFTWARE = "software"


class Direction(str, Enum):
    """Define constraints for direction"""

    SUPPORTS = "supports"
    NEUTRAL = "neutral"
    DISPUTES = "disputes"


#########################################
# GKS Common Abstract Utility Classes
# These do not inherit from Entity and are not typed explicitly
#########################################

class Code(RootModel):
    """Indicates that the value is taken from a set of controlled strings defined
    elsewhere. Technically, a code is restricted to a string which has at least one
    character and no leading or  trailing whitespace, and where there is no whitespace
    other than single spaces in the contents."""

    root: Annotated[str, StringConstraints(pattern=r'\S+( \S+)*')] = Field(
        ...,
        json_schema_extra={
            'description': 'Indicates that the value is taken from a set of controlled strings defined elsewhere. Technically, a code is restricted to a string which has at least one character and no leading or  trailing whitespace, and where there is no whitespace other than single spaces in the contents.',
            'example': 'ENSG00000139618',
        }
    )


class IRI(RootModel):
    """An IRI Reference (either an IRI or a relative-reference), according to `RFC3986
    section 4.1 <https://datatracker.ietf.org/doc/html/rfc3986#section-4.1>` and
    `RFC3987 section 2.1 <https://datatracker.ietf.org/doc/html/rfc3987#section-2.1>`.
    MAY be a JSON Pointer as an IRI fragment, as described by `RFC6901 section 6
    <https://datatracker.ietf.org/doc/html/rfc6901#section-6>`.
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
        json_schema_extra={'description': 'An IRI Reference (either an IRI or a relative-reference), according to `RFC3986 section 4.1  <https://datatracker.ietf.org/doc/html/rfc3986#section-4.1>` and `RFC3987 section 2.1 <https://datatracker.ietf.org/doc/html/rfc3987#section-2.1>`. MAY be a JSON Pointer as an IRI fragment, as  described by `RFC6901 section 6 <https://datatracker.ietf.org/doc/html/rfc6901#section-6>`.',
        }
    )



class Coding(BaseModel):
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
    code: Code = Field(
        ...,
        description="A symbol uniquely identifying the concept, as in a syntax defined by the code system. CURIE format is preferred where possible (e.g. 'SO:0000704' is the CURIE form of the Sequence Ontology code for 'gene')."
    )


class ConceptMapping(BaseModel):
    """A mapping to a concept in a terminology or code system."""

    model_config = ConfigDict(use_enum_values=True)

    coding: Coding = Field(..., description="A structured representation of a code for a defined concept in a terminology or code system.")
    relation: Relation = Field(..., description="A mapping relation between concepts as defined by the Simple Knowledge Organization System (SKOS).")


class Extension(BaseModel):
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


#########################################
# GKS Common Abstract Entity Class Definitions
#########################################

class Entity(BaseModel):
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
        description='A primary label for the entity.'
    )
    description: Optional[str] = Field(
        None,
        description='A free-text description of the Entity.'
    )
    alternativeLabels: Optional[List[str]] = Field(None, description="Alternative name(s) for the Entity.")
    extensions: Optional[List[Extension]] = Field(None, description="A list of extensions to the Entity, that allow for capture of information not directly supported by elements defined in the model.")


class DomainEntity(Entity):
    """An Entity that is specific to a particular biomedical domain such as disease,
    therapeutics, or genes. Domain Entities are considered as 'concept-level' entities,
    as opposed to particular instances. e.g. 'Lung Cancer', not 'patient123's lung
    cancer'. Or 'Erlotinib', not the particular doses given to a patient on a specific
    occasion.

    Abstract base class to be extended by other classes. Do NOT instantiate directly.
    """

    mappings: Optional[List[ConceptMapping]] = Field(None, description="A list of mappings to concepts in terminologies or code systems. Each mapping should include a coding and a relation.")


class Agent(Entity):
    """An autonomous actor (person, organization, or software agent) that bears some
    form of responsibility for an activity taking place, for the existence of an entity,
    or for another agent's activity.
    """

    type: Literal["Agent"] = Field(CoreImType.AGENT.value, description=f"MUST be '{CoreImType.AGENT.value}'.")
    name: Optional[str] = Field(None, description="The descriptive name of the agent.")
    subtype: Optional[AgentSubtype] = Field(None, description="A specific type of agent the Agent object represents.")


class ActivityBase(Entity):
    """Internal base class that holds shared fields for Activity model.

    Abstract base class to be extended by other classes. Do NOT instantiate directly.
    """

    subtype: Optional[Coding] = Field(None, description="A specific type of activity the Activity instance represents.")
    date: Optional[str] = Field(None, description="The date that the Activity was completed.")
    specifiedBy: Optional[List[Method]] = Field(None, description="A method that was followed in performing an Activity, that describes how it was executed.")

    @field_validator("date")
    @classmethod
    def date_format(cls, v: Optional[str]) -> Optional[str]:
        """Check that date is YYYY-MM-DD format"""
        if v:
            valid_format = "%Y-%m-%d"

            try:
                datetime.datetime.strptime(v, valid_format).replace(
                    tzinfo=datetime.timezone.utc
                ).strftime(valid_format)
            except ValueError:
                logging.warning("`date` SHOULD be formatted as a date string in ISO format 'YYYY-MM-DD'")
        return v

class Activity(ActivityBase):
    """An action or set of actions performed by an agent, that occurs over a period of
    time. Activities may use, generate, modify, move, or destroy one or more entities.
    """

    performedBy: Optional[List[Agent] ]= Field(None, description="An Agent who contributed to executing the Activity.")


class Contribution(ActivityBase):
    """An action taken by an agent in contributing to the creation, modification,
    assessment, or deprecation of a particular entity (e.g. a Statement, EvidenceLine,
    DataSet, Publication, etc.)
    """

    type: Literal["Contribution"] = Field(CoreImType.CONTRIBUTION.value, description=f"MUST be {CoreImType.CONTRIBUTION.value}.")
    contributor: Optional[List[Agent]] = Field(None, description="The agent that made the contribution.", min_length=1, max_length=1)
    activityType: Optional[Coding] = Field(None, description="The specific type of activity performed or role played by an agent in making the contribution (e.g. for a publication, agents may contribute as a primary author, editor, figure designer, data generator, etc. . Values of this property may be framed as activities or as contribution roles (e.g. using terms from the Contribution Role Ontology (CRO)).")


class InformationEntityBase(Entity):
    """Internal base class that holds shared fields for InformationEntity model.

    Abstract base class to be extended by other classes. Do NOT instantiate directly.
    """

    type: Literal["InformationEntity"] = Field(CoreImType.INFORMATION_ENTITY.value, description=f"MUST be {CoreImType.INFORMATION_ENTITY.value}.")
    specifiedBy: Optional[Union[Method, IRI]] = Field(None, description="A specification that describes all or part of the process that led to creation of the Information Entity ")
    contributions: Optional[List[Contribution] ]= Field(None, description="Specific actions taken by an Agent toward the creation, modification, validation, or deprecation of an Information Entity.")
    reportedIn: Optional[List[Union[Document, IRI]]] = Field(None, description="A document in which the the Information Entity is reported.")
    dateAuthored: Optional[str] = Field(None, description="Indicates when the information content expressed in the Information Entity was generated.")
    recordMetadata: Optional[RecordMetadata] = Field(None, description="Provenance metadata about a specific concrete record of information as encoded/serialized in a particular data set or object (as opposed to provenance about the abstract information content the encoding carries).")


class InformationEntity(InformationEntityBase):
    """An abstract (non-physical) entity that is about something - representing the
    underlying 'information content' conveyed by physical or digital information
    artifacts like books, web pages, data tables, or photographs.
    """

    derivedFrom: Optional[List[InformationEntity]] = Field(None, description="Another Information Entity from which this Information Entity is derived, in whole or in part.")

class Document(InformationEntity):
    """A collection of information, usually in a text-based or graphic human-readable
    form, intended to be read and understood together as a whole.
    """

    type: Literal["Document"] = Field(CoreImType.DOCUMENT.value, description=f"Must be '{CoreImType.DOCUMENT.value}'.")
    subtype: Optional[Coding] = Field(
        None, description="A specific type of document that a Document instance represents (e.g.  'publication', 'patent', 'pathology report')"
    )
    title: Optional[str] = Field(None, description="The official title given to the document by its authors.")
    urls: Optional[List[Annotated[str, StringConstraints(pattern=r"^(https?|s?ftp)://")]]] = Field(
        None, description="One or more URLs from which the content of the Document can be retrieved."
    )
    doi: Optional[Annotated[str, StringConstraints(pattern=r"^10\.(\d+)(\.\d+)*\/[\w\-\.]+")]] = Field(
        None,
        description="A `Digital Object Identifier <https://www.doi.org/the-identifier/what-is-a-doi/>_` for the document.",
    )
    pmid: Optional[int] = Field(
        None,
        description="A `PubMed unique identifier <https://en.wikipedia.org/wiki/PubMed#PubMed_identifier>`_.",
    )


class Method(InformationEntity):
    """A set of instructions that specify how to achieve some objective."""

    type: Literal["Method"] = Field(CoreImType.METHOD.value, description=f"MUST be '{CoreImType.METHOD.value}'.")
    subtype: Optional[Coding] = Field(
        None,
        description="A specific type of method that a Method instance represents (e.g. 'Variant Interpretation Guideline', or 'Experimental Protocol').",
    )
    license: Optional[str] = Field(None, description="A specific license that dictates legal permissions for how a method can be used (by whom, where, for what purposes, with what additional requirements, etc.).")


class RecordMetadata(BaseModel):
    """A reusable structure that encapsulates provenance metadata about a serialized
    data record or object in a particular dataset (as opposed to provenance about the
    real world entity this
    record or object represents).
    """

    recordIdentifier: Optional[str] = Field(None, description="The identifier of the data record or object described in this RecordMetadata object.")
    recordVersion: Optional[str] = Field(None, description="The version number of the record-level artifact the object describes.")
    derivedFrom: Optional[str] = Field(None, description="Another data record from which the record described here was derived, through a data ingest and/or transformation process. Value should be a string representing the identifier of the source record.")
    dateRecordCreated: Optional[str] = Field(None, description="The date the record was initially created.")
    contributions: Optional[List[Contribution]] = Field(None, description="Describes specific contributions made by an human or software agent to the creation, modification, or administrative management of a data record or object.")


class DataSet(InformationEntity):
    """A collection of related data items or records that are organized together in a
    common format or structure, to enable their computational manipulation as a unit.
    """

    type: Literal["DataSet"] = Field(CoreImType.DATA_SET.value, description=f"MUST be '{CoreImType.DATA_SET.value}'")
    subtype: Optional[Coding] = Field(None, description="A specific type of data set the DataSet instance represents (e.g. a 'clinical data set', a 'sequencing data set', a 'gene expression data set', a 'genome annotation data set')")
    releaseDate: Optional[str] = Field(None, description="Indicates when a version of a Data Set was formally released.")
    version: Optional[str] = Field(None, description="The version of the Data Set, as assigned by its creator.")
    license: Optional[str] = Field(None, description="A specific license that dictates legal permissions for how a data set can be used (by whom, where, for what purposes, with what additional requirements, etc.)")

class EvidenceLine(InformationEntity):
    """An independent, evidence-based argument that may support or refute the validity
    of a specific proposition. The strength and direction of this argument is based on
    an interpretation of one or more pieces of information as evidence for or against
    the target proposition.
    """

    type: Literal["EvidenceLine"] = Field(CoreImType.EVIDENCE_LINE.value, description=f"MUST be '{CoreImType.EVIDENCE_LINE.value}'")
    hasEvidenceItems: Optional[List[InformationEntity]] = Field(None, description="An individual piece of information that was evaluated as evidence in building the argument represented by an Evidence Line.")
    directionOfEvidenceProvided: Optional[Direction] = Field(None, description="The direction of support that the Evidence Line is determined to provide toward its target Proposition (supports, disputes, neutral)")
    strengthOfEvidenceProvided: Optional[Union[Coding, IRI]] = Field(None, description="The strength of support that an Evidence Line is determined to provide for or against its target Proposition, evaluated relative to the direction indicated by the directionOfEvidenceProvided value.")
    scoreOfEvidenceProvided: Optional[float] = Field(None, description="A quantitative score indicating the strength of support that an Evidence Line is determined to provide for or against its target Proposition, evaluated relative to the direction indicated by the directionOfEvidenceProvided value.")


class StatementBase(InformationEntity):
    """Internal base class that holds shared fields for Statement model.

    Abstract base class to be extended by other classes. Do NOT instantiate directly.
    """

    direction: Optional[Direction] = Field(None, description="A term indicating whether the Statement supports, disputes, or remains neutral w.r.t. the validity of the Proposition it evaluates.")
    strength: Optional[Union[Coding, IRI]]= Field(None, description="A term used to report the strength of a Proposition's assessment in the direction indicated (i.e. how strongly supported or disputed the Proposition is believed to be).  Implementers may choose to frame a strength assessment in terms of how *confident* an agent is that the Proposition is true or false, or in terms of the *strength of all evidence* they believe supports or disputes it.")
    score: Optional[float] = Field(None, description="A quantitative score that indicates the strength of a Proposition's assessment in the direction indicated (i.e. how strongly supported or disputed the Proposition is believed to be).  Depending on its implementation, a score may reflect how *confident* that agent is that the Proposition is true or false, or the *strength of evidence* they believe supports or disputes it.")
    statementText: Optional[str] = Field(None, description="A natural-language expression of what a Statement asserts to be true.")
    classification: Optional[Union[Coding, IRI]] = Field(None, description="A single term or phrase summarizing the outcome of direction and strength assessments of a Statement's proposition, in terms of a classification of its subject.")
    hasEvidenceLines: Optional[List[EvidenceLine]] = Field(None, description="An evidence-based argument that supports or disputes the validity of the proposition that a Statement assesses or puts forth as true. The strength and direction of this argument (whether it supports or disputes the proposition, and how strongly) is based on an interpretation of one or more pieces of information as evidence (i.e. 'Evidence Items).")



class Statement(StatementBase):
    """A claim of purported truth as made by a particular agent, on a particular
    occasion. Statements may be used to simply put forth a possible fact (i.e. a
    'proposition') as true, or to provide a more nuanced assessment of the level of
    confidence or evidence supporting a particular proposition.
    """

    subject: Dict = Field(..., description="The Entity about which the Statement is made.")
    predicate: str = Field(..., description="The relationship declared to hold between the subject and the object of the Statement.")
    object: Dict = Field(..., description="An Entity or concept that is related to the subject of a Statement via its predicate.")


class StudyGroup(Entity):
    """A collection of individuals or specimens from the same taxonomic class, selected
    for analysis in a scientific study based on their exhibiting one or more common
    characteristics  (e.g. species, race, age, gender, disease state, income). May be
    referred to as a 'cohort' or 'population' in specific research settings.
    """

    type: Literal["StudyGroup"] = Field(CoreImType.STUDY_GROUP.value, description=f"Must be '{CoreImType.STUDY_GROUP.value}'")
    memberCount: Optional[int] = Field(None, description="The total number of individual members in the StudyGroup.")
    isSubsetOf: Optional[List[StudyGroup] ]= Field(None, description="A larger StudyGroup of which this StudyGroup represents a subset.")
    characteristics: Optional[List[Characteristic]] = Field(None, description="A feature or role shared by all members of the StudyGroup, representing a criterion for membership in the group.")


class Characteristic(BaseModel):
    """An object holding a name-value pair used to describe a trait or role of an
    individual member of a StudyGroup.
    """

    name: str = Field(..., description="The type of the trait  or role described by the trait (e.g. 'ethnicity', 'sex', 'age', 'disease status').")
    value: str = Field(..., description="The specific value(s) that the indicated traitor role holds in all population members (e.g. 'east asian', 'female', 'adolescent', 'cancer').")
    valueOperator: Optional[bool] = Field(None, description="An operation that defines how to logically interpret a set of more than one Characteristic values ('AND', 'OR', 'NOT')")


class StudyResultBase(InformationEntityBase):
    """Internal base class that holds shared fields for StudyResult model.

    Abstract base class to be extended by other classes. Do NOT instantiate directly.
    """

    sourceDataSet: Optional[List[DataSet]] = Field(None, description="A larger DataSet from which the content of the StudyResult was derived.", max_length=1)
    ancillaryResults: Optional[Dict] = None
    qualityMeasures: Optional[Dict] = None


class StudyResult(InformationEntityBase):
    """A collection of data items from a single study that pertain to a particular
    subject or experimental unit in the study, along with optional provenance
    information describing how these data items were generated.

    Abstract base class to be extended by other classes. Do NOT instantiate directly.
    """

    focus: Optional[Union[DomainEntity, Coding, IRI]] = Field(None, description="The specific subject or experimental unit in a Study that data in the StudyResult object is about - e.g. a particular variant in a population allele frequency dataset like ExAC or gnomAD.")
    sourceDataSet: Optional[List[DataSet]] = Field(None, description="A larger DataSet from which the content of the StudyResult was derived.", max_length=1)
    componentResult: Optional[List[StudyResult]] = Field(None, description="Another StudyResult comprised of data items about the same focus as its parent Result, but based on a more narrowly scoped analysis of the foundational data (e.g. an analysis based on data about a subset of the parent Results full study population) .")
    studyGroup: Optional[StudyGroup] = Field(None, description="A description of a specific group or population of subjects interrogated in the ResearchStudy that produced the data captured in the StudyResult.")
    ancillaryResults: Optional[Dict] = None
    qualityMeasures: Optional[Dict] = None
