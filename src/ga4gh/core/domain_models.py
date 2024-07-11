"""GKS Common Library Domain Entity models

**This module should not be imported directly.**

Instead, users should use one of the following:

  * `from ga4gh.core import domain_models`, and refer to models with the
    abbreviated name, e.g., `domain_models.Gene` (recommended)

  * `import ga4gh.core`, and refer to models using the fully-qualified
    module name, e.g., `ga4gh.core.domain_models.Gene`
"""
from typing import Literal, Union, List
from enum import Enum

from pydantic import Field, RootModel

from ga4gh.core.entity_models import _DomainEntity


class CommonDomainType(str, Enum):
    """Define GKS Common Domain Entity types"""

    PHENOTYPE = "Phenotype"
    DISEASE = "Disease"
    TRAIT_SET = "TraitSet"
    TR_ACTION = "TherapeuticAction"
    TR_AGENT = "TherapeuticAgent"
    TR_SUB = "TherapeuticSubstituteGroup"
    TR_COMB = "CombinationTherapy"
    GENE = "Gene"

class Phenotype(_DomainEntity):
    """An observable characteristic or trait of an organism."""

    type: Literal[CommonDomainType.PHENOTYPE] = Field(
        CommonDomainType.PHENOTYPE,
        description=f'MUST be "{CommonDomainType.PHENOTYPE.value}".'
    )


class Disease(_DomainEntity):
    """A particular abnormal condition that negatively affects the structure or function
    of all or part of an organism and is not immediately due to any external injury.
    """

    type: Literal[CommonDomainType.DISEASE] = Field(
        CommonDomainType.DISEASE,
        description=f'MUST be "{CommonDomainType.DISEASE.value}".'
    )


class TraitSet(_DomainEntity):
    """A set of phenotype and/or disease concepts that together constitute a condition."""

    type: Literal[CommonDomainType.TRAIT_SET] = Field(
        CommonDomainType.TRAIT_SET,
        description=f'MUST be "{CommonDomainType.TRAIT_SET.value}".'
    )
    traits: List[Union[Disease, Phenotype]] = Field(
        ...,
        min_length=2
    )


class Condition(RootModel):
    """A disease or other medical disorder."""

    root: Union[TraitSet, Disease, Phenotype] = Field(
        ...,
        json_schema_extra={'description': 'A disease or other medical disorder.'},
        discriminator='type',
    )


class TherapeuticAction(_DomainEntity):
    """A therapeutic action taken that is intended to alter or stop a pathologic process."""

    type: Literal[CommonDomainType.TR_ACTION] = Field(
        CommonDomainType.TR_ACTION,
        description=f'MUST be "{CommonDomainType.TR_ACTION.value}".'
    )


class TherapeuticAgent(_DomainEntity):
    """An administered therapeutic agent that is intended to alter or stop a pathologic process."""

    type: Literal[CommonDomainType.TR_AGENT] = Field(
        CommonDomainType.TR_AGENT,
        description=f'MUST be "{CommonDomainType.TR_AGENT.value}".'
    )


class TherapeuticSubstituteGroup(_DomainEntity):
    """A group of therapeutic procedures that may be treated as substitutes for one another."""

    type: Literal[CommonDomainType.TR_SUB] = Field(
        CommonDomainType.TR_SUB,
        description=f'MUST be "{CommonDomainType.TR_SUB.value}".'
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

    type: Literal[CommonDomainType.TR_COMB] = Field(
        CommonDomainType.TR_COMB,
        description=f'MUST be "{CommonDomainType.TR_COMB.value}".'
    )
    components: List[Union[TherapeuticSubstituteGroup, TherapeuticAction, TherapeuticAgent]] = Field(
        ...,
        description='The individual therapeutic procedure components that constitute the combination therapy.',
        min_length=2
    )


class TherapeuticProcedure(RootModel):
    """An action or administration of therapeutic agents to produce an effect  that is
    intended to alter or stop a pathologic process.
    """

    root: Union[CombinationTherapy, TherapeuticSubstituteGroup, TherapeuticAction, TherapeuticAgent] = Field(
        ...,
        json_schema_extra={'description': 'An action or administration of therapeutic agents to produce an effect that is intended to alter or stop a pathologic process.'},
        discriminator='type',
    )


class Gene(_DomainEntity):
    """A basic physical and functional unit of heredity."""

    type: Literal[CommonDomainType.GENE] = Field(
        CommonDomainType.GENE,
        description=f'MUST be "{CommonDomainType.GENE.value}".'
    )
