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

The test suite includes VRS compliance tests from the submodule. Set the SeqRepo
path when running or recording tests, as CI does:

```shell
SEQREPO_ROOT_DIR=tests/data/seqrepo/latest python3 -m pytest --vcr-record=none
```

### Recorded test data

`pytest` and `make test` default to offline replay of recorded data, so they do not
need a SeqRepo REST service, a UTA database, or network access. Set
`SEQREPO_ROOT_DIR=tests/data/seqrepo/latest` (as CI does) when running or recording tests.

- **VCR cassettes** (`tests/**/cassettes/*.yaml`) replay vrs-python's HTTP requests to
  the [SeqRepo REST service](https://github.com/biocommons/seqrepo-rest-service).
  The default `--vcr-record=none` matches CI. In that mode (without `--disable-vcr`),
  opening a Python socket connection or connecting through psycopg2 fails the test
  immediately, including HTTP requests from tests not marked `@pytest.mark.vcr`
  (`block_live_connections` in `tests/conftest.py`). Connections opened by other C
  libraries are not intercepted.
- **The hgvs cache** (`tests/data/hgvs_cache.json`) replays the hgvs library's UTA
  queries and sequence fetches. Lookups that found no data are recorded too, and
  replay as the same `HGVSDataNotAvailableError`. A lookup that was never recorded
  fails immediately, even when production code would normally catch the HGVS error.
  `make record-hgvs-cache`
  uses `--record-hgvs-cache=PATH --disable-vcr` to record live data
  to a temporary pickle, then converts it to JSON if tests pass. hgvs itself reads
  and writes pickles; `tests/hgvs_cache_json.py` converts between the two formats.

To re-record, start the external services (`docker compose up` starts
seqrepo-rest-service on port 5000 and UTA on port 5432). With an existing local SeqRepo
snapshot, you can instead run just the REST service with it mounted:
`docker run -p 5000:5000 -v /path/to/seqrepo:/usr/local/share/seqrepo:ro biocommons/seqrepo-rest-service:0.2.2 seqrepo-rest-service /usr/local/share/seqrepo/2024-12-20`.
When a change needs new data in both, re-record the hgvs cache first, since cassette
recording reads hgvs data from it:

```shell
# hgvs cache: needs UTA and seqrepo-rest-service; does not touch cassettes
export UTA_DB_URL=postgresql://anonymous@localhost:5432/uta/uta_20241220
make record-hgvs-cache

# VCR cassettes: needs seqrepo-rest-service
make clean-cassettes
pytest --vcr-record=once
```

`--disable-vcr` allows live connections and bypasses HTTP cassette replay and writing,
even with the default `--vcr-record=none`. HGVS still replays its cache unless
`--record-hgvs-cache=PATH` is also supplied. The recording target uses both flags,
so it records HGVS data without changing HTTP cassettes.

`SEQREPO_REST_URL` (default `http://localhost:5000/seqrepo`) sets the REST service used
by both.

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
