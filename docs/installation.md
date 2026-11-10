# Installation

VRS-Python requires **Python 3.11 or later**. Install it in a virtual environment:

```shell
python -m venv venv
source venv/bin/activate
python -m pip install ga4gh.vrs
```

The base package supports model validation, serialization, and identifier computation.

## Translation, local sequence access, and VCF annotation

Install `extras` for [SeqRepo](https://github.com/biocommons/seqrepo),
[UTA](https://github.com/biocommons/uta),
[HGVS support](https://github.com/biocommons/hgvs), and VCF annotation.

### Prerequisites

Install the prerequisites for your platform:

#### macOS

With [Homebrew](https://docs.brew.sh/Installation) installed, run:

```shell
brew update
brew install libpq
brew install python3
brew install postgresql@14
```

#### Ubuntu

```shell
sudo apt install gcc libpq-dev python3-dev
```

### Install the extras

In your virtual environment, run:

```shell
python -m pip install 'ga4gh.vrs[extras]'
```

Sequence and transcript datasets require separate setup:

- Set up SeqRepo for allele normalization, SPDI/gnomAD translation, and VCF annotation.
- Set up UTA for HGVS workflows that need transcript metadata or alignments.
- Model construction, serialization, and identifier computation can run without
  either service; reference sequence validation requires SeqRepo.

Follow [Set up SeqRepo and UTA](data-services.md) for configuration instructions.

## Notebooks

The `notebooks` group includes dependencies for running the
[example notebooks](notebooks/index.md). Include `extras` for the translation and
sequence-access examples:

```shell
python -m pip install 'ga4gh.vrs[extras,notebooks]'
```
