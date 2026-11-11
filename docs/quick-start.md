# Quick start

Build a VRS allele, compute its identifier, and translate between HGVS and VRS.

- **Build and identify:** uses SeqRepo to resolve a reference sequence accession.
- **Translate HGVS:** also configures UTA for transcript metadata and alignments.

Identifier computation itself does not require an external service.

## Build and identify an allele

Before running this example:

1. Install [extras](installation.md#install-the-extras).
2. Start [SeqRepo](data-services.md#docker-compose).

Create an allele for `NC_000019.10:g.44908684T>C`:

```python
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import create_dataproxy

dp = create_dataproxy("seqrepo+http://localhost:5000/seqrepo")

# Resolve the RefSeq accession to a RefGet accession.
refget_accession = dp.derive_refget_accession("NC_000019.10")
if refget_accession is None:
    raise ValueError("NC_000019.10 was not found in the SeqRepo dataset")

allele = models.Allele(
    # Reference sequence and interval to replace.
    location=models.SequenceLocation(
        sequenceReference=models.SequenceReference(
            refgetAccession=refget_accession
        ),
        # Zero-based, inter-residue coordinates spanning residue 44908684.
        start=44908683,
        end=44908684,
    ),
    state=models.LiteralSequenceExpression(sequence="C"),  # Replacement sequence.
)
```

Compute the identifier and print the allele as JSON:

```python
identifier = ga4gh_identify(allele)
print(identifier)
print(allele.model_dump_json(exclude_none=True, indent=2))
```

The identifier is:

```text
ga4gh:VA.NsWfIPdbQ9LWo4Z5f9lIhMW87kZCtrEc
```

In the JSON, `id` is the GA4GH identifier, `digest` is its digest component,
and `type` names the VRS object type:

```json
{
  "id": "ga4gh:VA.NsWfIPdbQ9LWo4Z5f9lIhMW87kZCtrEc",
  "type": "Allele",
  "digest": "NsWfIPdbQ9LWo4Z5f9lIhMW87kZCtrEc",
  "location": {
    "type": "SequenceLocation",
    "digest": "bqt0x0LbTvEpKXw2x7cUfww7mI_Q_pZJ",
    "sequenceReference": {
      "type": "SequenceReference",
      "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    },
    "start": 44908683,
    "end": 44908684
  },
  "state": {
    "type": "LiteralSequenceExpression",
    "sequence": "C"
  }
}
```

## Translate HGVS to VRS

Start [UTA](data-services.md#docker-compose) before running this example.

Configure UTA and create a translator with the data proxy from above:

```python
import os

os.environ["UTA_DB_URL"] = "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20241220"

from ga4gh.vrs.extras.translator import AlleleTranslator

translator = AlleleTranslator(dp)
```

Translate the same HGVS expression. By default, translation
[normalizes the allele](api/normalization.md) and computes its identifier:

```python
allele = translator.to_vrs("NC_000019.10:g.44908684T>C", fmt="hgvs")
print(allele.model_dump_json(exclude_none=True, indent=2))
```

The allele and its identifier match the first example. Translation also assigns
an `id` to the location:

```json
{
  "id": "ga4gh:VA.NsWfIPdbQ9LWo4Z5f9lIhMW87kZCtrEc",
  "type": "Allele",
  "digest": "NsWfIPdbQ9LWo4Z5f9lIhMW87kZCtrEc",
  "location": {
    "id": "ga4gh:SL.bqt0x0LbTvEpKXw2x7cUfww7mI_Q_pZJ",
    "type": "SequenceLocation",
    "digest": "bqt0x0LbTvEpKXw2x7cUfww7mI_Q_pZJ",
    "sequenceReference": {
      "type": "SequenceReference",
      "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    },
    "start": 44908683,
    "end": 44908684
  },
  "state": {
    "type": "LiteralSequenceExpression",
    "sequence": "C"
  }
}
```

## Translate back to HGVS

`from_vrs` returns a list because a sequence can have multiple accessions.
By default, it uses RefSeq accessions in your SeqRepo dataset:

```python
hgvs_expressions = translator.from_vrs(allele, fmt="hgvs")
print(hgvs_expressions)
```

The result includes the original expression; other entries depend on your dataset:

```text
['NC_000019.10:g.44908684T>C']
```

See the [translation API](api/translation-and-hgvs.md) for options or
[explore the notebooks](notebooks/index.md) for more examples.
