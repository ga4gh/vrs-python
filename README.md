# VRS-Python

**VRS-Python** provides Python language support and a reference implementation for the
[GA4GH Variation Representation Specification(VRS)](https://github.com/ga4gh/vrs).

See the [**documentation**](https://ga4gh.github.io/vrs-python/latest/) for installation, usage, and API reference.

## Information

[![license](https://img.shields.io/badge/license-Apache-green)](https://github.com/ga4gh/vrs-python/blob/main/LICENSE) [![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ga4gh/vrs-python/main?labpath=notebooks)

## Releases

[![gitHub tag](https://img.shields.io/github/v/tag/ga4gh/vrs-python.svg)](https://github.com/ga4gh/vrs-python/releases) [![pypi](https://img.shields.io/pypi/v/ga4gh.vrs.svg)](https://pypi.org/project/ga4gh.vrs/) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14013256.svg)](https://doi.org/10.5281/zenodo.14013256)

## Development

 [![action status](https://github.com/ga4gh/vrs-python/actions/workflows/python-cqa.yml/badge.svg)](https://github.com/ga4gh/vrs-python/actions/workflows/python-cqa.yml) [![issues](https://img.shields.io/github/issues-raw/ga4gh/vrs-python.svg)](https://github.com/ga4gh/vrs-python/issues)
[![GitHub Open Pull Requests](https://img.shields.io/github/issues-pr/ga4gh/vrs-python.svg)](https://github.com/ga4gh/vrs-python/pull/) [![GitHub Contributors](https://img.shields.io/github/contributors/ga4gh/vrs-python.svg)](https://github.com/ga4gh/vrs-python/graphs/contributors/) [![GitHub stars](https://img.shields.io/github/stars/ga4gh/vrs-python.svg?style=social&label=Stars)](https://github.com/ga4gh/vrs-python/stargazers) [![GitHub forks](https://img.shields.io/github/forks/ga4gh/vrs-python.svg?style=social&label=Forks)](https://github.com/ga4gh/vrs-python/forks)

## Features

- Pydantic implementation of GKM-Core models and VRS models
- Algorithm for generating consistent, globally unique identifiers for variation without a central authority
- Algorithm for performing fully justified allele normalization
- Translating from and to other variant formats
- Annotate VCFs with VRS
- Convert GA4GH objects between inlined and referenced forms

## Security Note (from the GA4GH Security Team)

A stand-alone security review has been performed on the specification itself.
This implementation is offered as-is, and without any security guarantees. It
will need an independent security review before it can be considered ready for
use in security-critical applications. If you integrate this code into your
application it is AT YOUR OWN RISK AND RESPONSIBILITY to arrange for a security
audit.
