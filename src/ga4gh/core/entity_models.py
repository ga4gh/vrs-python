"""GKS Common Library Entity models

**This module should not be imported directly.**

Instead, users should use one of the following:

  * `from ga4gh.core import entity_models`, and refer to models with the
    abbreviated name, e.g., `entity_models.Coding` (recommended)

  * `import ga4gh.core`, and refer to models using the fully-qualified
    module name, e.g., `ga4gh.core.entity_models.Coding`
"""
from abc import ABC
from typing import Any, Dict, Annotated, Optional, Union, List
from enum import Enum

from pydantic import BaseModel, Field, RootModel, StringConstraints, ConfigDict

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
    version: Optional[str] = Field(
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
        None, description='The value of the Extension - can be any primitive or structured object'
    )
    description: Optional[str] = Field(None, description="A description of the meaning or utility of the Extension, to explain the type of information it is meant to hold.")


class Expression(BaseModel):
    """Representation of a variation by a specified nomenclature or syntax for a
    Variation object. Common examples of expressions for the description of molecular
    variation include the HGVS and ISCN nomenclatures.
    """

    model_config = ConfigDict(use_enum_values=True)

    syntax: Syntax = Field(..., description="The syntax used to describe the variation. The value should be one of the supported syntaxes.")
    value: str = Field(..., description="The expression of the variation in the specified syntax. The value should be a valid expression in the specified syntax.")
    syntax_version: Optional[str] = Field(None, description="The version of the syntax used to describe the variation. This is particularly important for HGVS expressions, as the syntax has evolved over time.")


#########################################
# GKS Common Abstract Entity Class Definitions
#########################################


class Entity(ABC, BaseModel):
    """Entity is the root class of the 'gks-common' core information model classes -
    those that have identifiers and other general metadata like labels, xrefs, urls,
    descriptions, etc. All common classes descend from and inherit its attributes.
    """

    id: Optional[str] = Field(
        None,
        description="The 'logical' identifier of the entity in the system of record, e.g. a UUID. This 'id' is unique within a given system. The identified entity may have a different 'id' in a different system, or may refer to an 'id' for the shared concept in another system (e.g. a CURIE)."
    )
    type: str
    label: Optional[str] = Field(
        None,
        description='A primary label for the entity.'
    )
    description: Optional[str] = Field(
        None,
        description='A free-text description of the entity.'
    )
    alternativeLabels: Optional[List[str]] = Field(None, description="Alternative name(s) for the Entity.")
    extensions: Optional[List[Extension]] = Field(None, description="A list of extensions to the entity. Extensions are not expected to be natively understood, but may be used for pre-negotiated exchange of message attributes between systems.")


class DomainEntity(Entity, ABC):
    """An Entity that is specific to a particular biomedical domain such as disease,
    therapeutics, or genes. Domain Entities are considered as 'concept-level' entities,
    as opposed to particular instances. e.g. 'Lung Cancer', not 'patient123's lung
    cancer'. Or 'Erlotinib', not the particular doses given to a patient on a specific
    occasion.
    """

    mappings: Optional[List[ConceptMapping]] = Field(None, description="A list of mappings to concepts in terminologies or code systems. Each mapping should include a coding and a relation.")
