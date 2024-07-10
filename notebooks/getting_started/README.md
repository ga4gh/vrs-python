# VRS Python Getting Started Notebook Series
This series of notebooks is intended to be the fastest way to becoming productive with the
[vrs-python](https://github.com/ga4gh/vrs-python) package. The
intent of each notebook in this series is to be an interactive introduction to functionally contained in the vrs-python package.

A beginning developer level of familiarity with python, jupyter notebooks is assumed in order to run this notebook series.
You should be familiar with installing packages and running commands in your execution environment.

### Pre-requisites
The following software packages must exist in your execution environment before running these notebooks:
* git
* python@3.12
* make

### Setup vrs-python
From a terminal window, run the following commands:
* git clone --recurse-submodules https://github.com/ga4gh/vrs-python
* cd vrs-python
* make nbready
* source venv/3.12/bin/activate
* cd notebooks/getting_started
* jupyter notebook notebook_name.ipynb

## Notebooks
### 1 Quick Start
The [Quick Start](1_Quick_Start.ipynb) notebook details how to get started by
setting up access to a SeqRepo *DataProxy* and introduces the user to using an *AlleleTranslator* to convert
the same allele to VRS form from both it's SPDI and HGVS nomenclature forms.
### 2 Exploring the SeqRepo Data Proxy
Sequence references are at the core of many of the operations for converting to and from VRS variant representations.
The [Exploring the SeqRepo Data Proxy](2_Exploring_the_SeqRepo_DataProxy.ipynb) notebook
has a number of useful utility methods for accessing information about sequence references.
### 3 Basic Models
In the [Basic Models](3_Basic_Models.ipynb) notebook, we explore building a VRS *Allele*
from its component parts. The notebook details how to add VRS identifiers to the identifiable components of the *Allele*.
### 4 Exploring the Allele Translator
The current implementation of vrs-python facilitates transformation of variants
in a number of different variant nomenclatures (SPDI, HGVS, gnomAD and Beacon) to VRS form. In the
[Exploring the Allele Translator](4_Exploring_the_AlleleTranslator.ipynb) notebook,
we show how to transform basic variants to VRS, and in some cases, back to the original nomenclature of the variant.
### 5 Exploring the CNV Translator
The final notebook of this series,
[Exploring the CNV Translator](5_Exploring_the_CnvTranslator.ipynb) details transformations
of various forms of copy number variation to their VRS representations.





