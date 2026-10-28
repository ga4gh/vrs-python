# DataProxy backends

VRS-Python accesses reference sequences and resolves cross-namespace sequence
identifiers (e.g. `refseq:NM_000551.3` → `ga4gh:SQ.…`) through a **DataProxy**.
A DataProxy combines two responsibilities:

- **Sequence access** — `get_sequence()`, `get_sequence_length()`
- **Alias / identifier resolution** — `get_aliases()`,
  `translate_sequence_identifier()`, `derive_refget_accession()`

You construct one with `create_dataproxy(uri)`, where the URI scheme selects the
backend:

```python
from ga4gh.vrs.dataproxy import create_dataproxy

dp = create_dataproxy("seqrepo+file:///usr/local/share/seqrepo/2024-12-20")
```

If `uri` is omitted, `create_dataproxy` falls back to the
`GA4GH_VRS_DATAPROXY_URI` environment variable.

## Available backends and URI schemes

| Scheme | Backend | Notes |
| --- | --- | --- |
| `seqrepo+file:///abs/path` | Local [SeqRepo](https://github.com/biocommons/biocommons.seqrepo) | Default backend; reads a SeqRepo directory. |
| `seqrepo+:relative/path` | Local SeqRepo | Same as above with a relative path (note the empty protocol). |
| `seqrepo+http://host:port/seqrepo` | [SeqRepo REST service](https://github.com/biocommons/seqrepo-rest-service) | Remote SeqRepo over HTTP. |
| `seqrepo+https://host/seqrepo` | SeqRepo REST service | TLS variant. |
| `refget+file:///abs/path` | [gtars](https://github.com/databio/gtars) RefgetStore | Local on-disk RefgetStore. Requires `pip install ga4gh.vrs[refget]`. |
| `refget+http://host[:port]/path` | gtars RefgetStore | Remote static HTTP(S) store; working set is cached locally. |
| `refget+https://host/path` | gtars RefgetStore | TLS variant. |

> **Relative paths:** `urlparse` treats the first segment after `file://` as a
> host, so `refget+file://my/store` resolves to the path `/store`. Use an
> absolute path (`refget+file:///abs/path`) or the empty-protocol form
> (`refget+:relative/path`).

The `refget` backend composes a `RefgetSequenceDataProxy` (sequence retrieval by
GA4GH digest) and a `RefgetAliasDataProxy` (alias resolution via the store's
alias index) behind a `CompositeDataProxy`: identifiers are tried against the
sequence store first, and on a miss are resolved to a `ga4gh:SQ.` digest via the
alias index and retried.

### Remote refget caching

Remote `refget+http(s)://` stores cache their index, alias sidecars, and fetched
sequence payloads on the local filesystem. The cache root resolves to, in order:

1. `GA4GH_VRS_REFGET_CACHE_DIR`, if set
2. `$XDG_CACHE_HOME/gtars/refget`, if `XDG_CACHE_HOME` is set
3. `~/.cache/gtars/refget`

## In-memory caching

Separately from the on-disk refget cache above, all backends memoize lookups
in-process with `functools.lru_cache`. Two environment variables tune the sizes:

| Variable | Default | Covers |
| --- | --- | --- |
| `GA4GH_VRS_DATAPROXY_CACHE_SIZE` | `8192` | Sequence lengths, digests, and alias lists, across all backends. |
| `GA4GH_VRS_REST_DATAPROXY_SEQUENCE_CACHE_SIZE` | `128` | Sequence payloads fetched by the SeqRepo REST backend. |

These are separate knobs because the cached values differ in size by orders of
magnitude. The first holds small fixed-size values — an integer length, a
digest, a list of a dozen short alias strings — so 8192 entries costs a few MB.
The second holds sequence text, which is bounded only by the caller's requested
span: the normalization path stores 1-2 base strings per entry, but
`get_sequence(identifier)` with no bounds stores an entire contig in one entry.

`GA4GH_VRS_REST_DATAPROXY_SEQUENCE_CACHE_SIZE` applies only to
`SeqRepoRESTDataProxy` — the only backend that caches sequence payloads at all.
Local SeqRepo and refget `get_sequence()` calls hit the store every time.

**This cache memoizes fetched regions, not sequences.** The key is the exact
`(identifier, start, end)` triple that was requested, so it is only useful when
the identical region is fetched repeatedly. In particular:

- An overlapping or adjacent span is a miss. Having `(chr1, 100, 200)` cached
  does nothing for a request for `(chr1, 150, 250)` or `(chr1, 200, 300)`.
- A cached full-sequence entry cannot serve a sub-range request. `get_sequence(id)`
  is stored under `(id, None, None)` and is never consulted for `(id, 0, 10)`.
- A linear scan over distinct positions gets no benefit at all — every request is
  a fresh key. Raising this value does not speed up streaming workloads such as
  annotating a VCF of distinct variants; it only helps genuinely repeated
  identical lookups (e.g. re-normalizing the same allele, or several alleles
  resolving the same span).
- Positional and keyword arguments hash to different keys, a `functools.lru_cache`
  property. `get_sequence(ac, 5, 9)` and `get_sequence(ac, start=5, end=9)` are
  two separate entries and two separate HTTP requests, even though they request
  the same bases.

So treat it as a de-duplication cache for repeated identical fetches, not as a
read-ahead or a local mirror of the sequence.

Accepted values for both:

| Value | Meaning |
| --- | --- |
| positive integer | Maximum number of cached entries |
| `0` | Disable caching (`functools` semantics) |
| `-1` or `none` | Unbounded cache |
| anything else | Raises `ValueError` at import |

> **Both are read once, at import time.** The `lru_cache` decorators are applied
> when the classes are defined, so setting these after `ga4gh.vrs.dataproxy` has
> been imported has no effect. Export them before starting the process.

## Known divergences and caveats (refget vs SeqRepo)

The refget and SeqRepo backends are functionally equivalent for the sequence and
alias data VRS-Python consumes, with two differences worth knowing about.

### Stop-codon (`*`) digest divergence on protein sequences

For protein sequences that contain embedded stop-codon `*` characters, a
RefgetStore can compute a **different `ga4gh:SQ.` digest than SeqRepo for the
same underlying sequence bytes**.

- **Cause:** SeqRepo runs sequences through
  `bioutils.sequences.normalize_sequence()`, which strips `*` before computing
  the `sha512t24u` digest. gtars computes the digest directly from the raw FASTA
  bytes, including any `*`.
- **Which is correct:** The RefgetStore digest. The GA4GH VRS specification does
  not prescribe `*` stripping, so hashing the raw bytes is spec-correct.
- **Scope:** Small. In a full Ensembl r113 + RefSeq load, 116 of 322,402
  comparable Ensembl aliases diverged (~0.04%); 115 were pure `*`-normalization
  artifacts (identical bytes, different digest) and 1 was genuine cross-release
  sequence drift. **NCBI RefSeq protein FASTAs contain no `*` and showed zero
  digest mismatches.** DNA and RNA sequences are never affected.

The practical implication: a VRS Allele on an affected Ensembl protein will have
a different sequence identifier depending on which backend produced it. If you
need identifiers that match an existing SeqRepo-derived corpus on these specific
accessions, stay on SeqRepo for protein annotation.

### Version-drift policy

SeqRepo is **accumulative**: it retains every version of every accession it has
ever loaded across historical imports. A RefgetStore typically pins **one
version** of each source (e.g. a single Ensembl release; the current NCBI RefSeq
transcript/protein shards, which publish only the latest version of each
accession). As a result, the refget backend may raise `KeyError` on older
transcript/protein versions that a long-lived SeqRepo instance still holds.
