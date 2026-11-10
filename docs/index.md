# VRS-Python

VRS-Python implements the
[GA4GH Variation Representation Specification (VRS)](https://vrs.ga4gh.org/)
in Python. It helps you represent variants consistently and compare records
across data sources using identifiers computed from their content.

Use VRS-Python to:

- Build and validate VRS and GKM-Core objects with Pydantic
  (see [VRS compatibility](compatibility.md)).
- Normalize alleles and compute consistent GA4GH identifiers.
- Translate supported HGVS, SPDI, and gnomAD expressions to VRS.
- Annotate VCF files with VRS allele identifiers.
- Store and exchange objects using inline or referenced components.

## Start using VRS-Python

<div class="grid cards" markdown>

- **[Install](installation.md)**

    Choose the dependencies for your workflow.

- **[Quick start](quick-start.md)**

    Build your first allele and compute its identifier.

- **[Explore with notebooks](notebooks/index.md)**

    Work through models, sequence access, and translation examples.

- **[API reference](api/index.md)**

    Find Python classes, functions, and their parameters.

- **[Set up SeqRepo and UTA](data-services.md)**

    Connect sequence and transcript data for translation and normalization.

- **[Contribute](development.md)**

    Improve the implementation, examples, and documentation.

</div>
