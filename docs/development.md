# Development and contributing

## Installing for development

Install Python 3.14, available as `python3.14`, before running `make devready`.

[Fork the GitHub repo](https://github.com/ga4gh/vrs-python/fork), then clone your
fork and initialize a development environment:

```shell
git clone --recurse-submodules git@github.com:YOUR_GITHUB_ID/vrs-python.git
cd vrs-python
make devready
source venv/3.14/bin/activate
```

`make devready` installs development dependencies and pre-commit hooks. For a
manually prepared environment, install the development dependencies:

```shell
python -m pip install -e '.[tests,lint,docs,extras,notebooks]'
```

If you already cloned the repo, but forgot to include `--recurse-submodules` you can run:

```shell
git submodule update --init --recursive
```

Keep the VRS submodule at the commit recorded by this repository when switching
branches. Do not update it incidentally while changing documentation or code.

## Set up SeqRepo and UTA

Depending on the code or notebooks you work on, you may need to configure
SeqRepo for reference sequences and UTA for transcript metadata and alignments.
See [Set up SeqRepo and UTA](data-services.md) for setup options and environment
variables.

## Install the pre-commit hooks

If you installed the development dependencies manually, install the hooks:

```shell
pre-commit install
```

This step is already included in `make devready`.

## Run the style checks

After making your changes, check linting and formatting:

```shell
python3 -m ruff check .
python3 -m ruff format --check .
```

## Run the tests

Add focused tests for changed behavior, then run the test suite:

```shell
python3 -m pytest
```

The test suite includes VRS compliance tests from the submodule. For direct pytest
runs, set the SeqRepo path and use recorded network responses, as CI does:

```shell
SEQREPO_ROOT_DIR=tests/data/seqrepo/latest python3 -m pytest --vcr-record=none
```

## Build and preview the documentation

Edit the documentation sources:

- Markdown pages in `docs`.
- API docstrings in `src`.
- Notebooks in `notebooks/getting_started`.

The build generates `docs-build` and `site`. Do not edit or commit these directories.

Build the documentation with strict checks before opening a pull request:

```shell
make docs
```

To preview the documentation locally:

```shell
make docs-serve
```

Open <http://127.0.0.1:8000>. After editing source pages or notebooks, stop the
server and rerun the preview command to refresh the staged documentation.

## Open a pull request

Run the checks relevant to your changes, review the generated documentation, and
open a pull request on GitHub against `main`. Include a short summary of the change
and any relevant testing or documentation notes.

- Use `Closes #123` to close an issue automatically when the pull request merges.
- Use `Refs #123` or link the issue in GitHub for related issues that remain open.
