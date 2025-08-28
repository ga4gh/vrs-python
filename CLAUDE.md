# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

VRS-Python is the reference implementation of the GA4GH Variation Representation Specification (VRS), providing Python language support for representing genomic variation data consistently and uniquely.

## Development Commands

### Environment Setup
- `make devready` - Create Python 3.12 venv, install dev dependencies, setup pre-commit hooks
- `make nbready` - Create Python 3.12 venv, install notebook dependencies for Jupyter usage
- `source venv/3.12/bin/activate` - Activate the virtual environment (required before development)

### Development Tasks
- `make test` - Run pytest test suite
- `make lint` - Static analysis with ruff (auto-fix enabled)
- `make format` - Code formatting with ruff
- `pytest tests/path/to/specific_test.py` - Run individual test files
- `pytest tests/path/to/specific_test.py::test_function_name` - Run specific test functions

### Docker Services
- `docker-compose up` - Start external dependencies (seqrepo-rest-service on port 5000, UTA database on port 5432)
- Required for full functionality of VRS translation and normalization features

### Package Installation
- `pip install -e .[dev,extras,notebooks]` - Install in development mode with all dependencies
- `pip install -e .[extras]` - Install with core extras (SeqRepo, HGVS tools, etc.)

## Code Architecture

### Core Components
- **`src/ga4gh/vrs/models.py`** - Pydantic models for VRS objects (Allele, Location, etc.)
- **`src/ga4gh/vrs/normalize.py`** - Allele normalization algorithms per VRS specification
- **`src/ga4gh/vrs/dataproxy.py`** - Abstract interface and SeqRepo implementation for sequence data access
- **`src/ga4gh/vrs/enderef.py`** - Converting between inlined and referenced VRS object forms

### Key Modules
- **`src/ga4gh/vrs/extras/translator.py`** - Translates between VRS and external formats (HGVS, SPDI, gnomAD, Beacon)
- **`src/ga4gh/vrs/extras/annotator/`** - VCF annotation tools with VRS identifiers
- **`src/ga4gh/vrs/utils/hgvs_tools.py`** - HGVS parsing and validation utilities
- **`src/ga4gh/core/`** - Core GA4GH models and identifier generation

### Testing Structure
- **`tests/`** - Main test directory with pytest configuration
- **`tests/validation/`** - VRS specification compliance tests
- **`tests/extras/`** - Tests for translator and annotator functionality
- **`tests/cassettes/`** - VCR.py cassettes for external API mocking

## Development Notes

### Dependencies
- Requires Python >= 3.10 (Python 3.12 required for development)
- External services: SeqRepo (sequence data), UTA (transcript alignments)
- Key libraries: Pydantic 2.x, bioutils, HGVS, requests

### Testing Environment
- Uses pytest with coverage reporting
- VCR.py for API response caching
- Test data includes local SeqRepo instance at `tests/data/seqrepo/`
- Set `SEQREPO_ROOT_DIR=tests/data/seqrepo/latest` for tests

### Pre-commit Configuration
- Ruff linting and formatting
- Automatic code quality checks before commits
- Install with `pre-commit install` (done automatically by `make develop`)
