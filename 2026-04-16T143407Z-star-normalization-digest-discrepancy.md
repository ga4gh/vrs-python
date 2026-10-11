# Stop Codon (`*`) Normalization Discrepancy in Protein Sequence Digests

**Date**: 2026-04-16
**Branch**: `kf/dataproxy-backend-refactor`

## 1. Summary

When loading Ensembl release 113 protein sequences into a gtars
RefgetStore, 115 ENSP accessions produce different `ga4gh:SQ.*`
identifiers than the same accessions in seqrepo 2024-12-20 — despite
the underlying sequence bytes being identical in both stores.

The root cause is that these 115 proteins contain internal `*` (stop
codon / Ter) characters in their FASTA sequences, and the two systems
handle `*` differently when computing the sequence digest:

- **SeqRepo** (via bioutils) strips `*` before hashing — conformant
  with the RefGet spec's allowed alphabet (A-Z only)
- **gtars RefgetStore** hashes the raw FASTA bytes including `*` —
  non-conformant with the RefGet spec for these sequences

This discrepancy was invisible until we loaded Ensembl protein data
because NCBI RefSeq protein FASTAs never contain `*` characters. Zero
of the 197,939 RefSeq protein sequences have `*`, so the normalization
difference was never triggered.

## 2. Which proteins are affected?

All 115 are Ensembl-only protein annotations — none has a RefSeq (NP_/
XP_) counterpart with the same sequence, even after stripping `*`.
These fall into a few gene families:

| Category | Count | Description |
|---|---|---|
| IG D-segments (IGHD*) | 26 | Immunoglobulin diversity segments — short (5-12 aa) peptides from D-region CDS that inherently hit stop codons in translation |
| KIR receptors (KIR2DS4) | 25 | Killer cell immunoglobulin receptors — highly polymorphic immune locus with frameshift/stop alleles |
| Olfactory receptors (OR*) | 22 | Many OR genes are pseudogenes with in-frame stops |
| CASP12 | 9 | Caspase-12 — pseudogene in most human populations (premature stop) |
| CYP2D7 | 8 | Cytochrome P450 pseudogene |
| Other | 25 | MUC6, PKD1L2, POTEA, IFNL4, GBA3, GSTT2, AKR7L, MROH5, etc. |

Ensembl annotates CDS regions for these genes even when translation
crosses stop codons. RefSeq's curation is more conservative: if a CDS
hits a stop, RefSeq either omits the protein record entirely, or picks
a different allele/isoform that avoids the stop. This is why the issue
surfaced only when we started loading Ensembl protein data.

### Stop codon positions

153 of the 158 `*`-containing Ensembl proteins have `*` at **internal**
positions (not at the end). Only 4 have terminal-only `*`, and 1 has
both. Examples from the FASTA:

```
>ENSP00000488168.1 gene_symbol:IGHD4-23 gene_biotype:IG_D_gene (non-functional)
*LRW*L              ← 6 aa, * at positions 1 and 5

>ENSP00000488711.1 gene_symbol:IGHD3-22 gene_biotype:IG_D_gene
VLL***WLLL          ← 10 aa, 3 consecutive stops at positions 4-6

>ENSP00000340011.8 gene_symbol:KIR2DS4 gene_biotype:protein_coding
...TEQ*TARILMNKTIRRCHTH  ← 297 aa, * at positions 239, 269, 280
```

The 115 mismatched proteins contain between 1 and 23 internal stop
codons:

| `*` count | Proteins |
|---|---|
| 1 | 54 |
| 2 | 19 |
| 3 | 17 |
| 4-6 | 9 |
| 7-12 | 10 |
| 13-23 | 6 |

## 3. How the digest is computed — seqrepo vs gtars

### SeqRepo (bioutils)

SeqRepo's `store` method imports `seq_seqhash` aliased as `sha512t24u`:

```python
# biocommons/seqrepo/seqrepo.py, line 11
from bioutils.digests import seq_seqhash as sha512t24u
```

When a sequence is stored (line ~148):

```python
seqhash = sha512t24u(seq)  # actually calls seq_seqhash(seq)
seq_id = seqhash
```

`seq_seqhash` normalizes before hashing:

```python
# bioutils/digests.py — seq_seqhash
def seq_seqhash(seq, normalize=True):
    seq = normalize_sequence(seq) if normalize else seq
    return str(vmc_digest(seq, digest_size=24))
```

And `normalize_sequence` strips `*`:

```python
# bioutils/sequences.py — normalize_sequence
def normalize_sequence(seq):
    """Removes whitespace and asterisks, and uppercases the string."""
    nseq = re.sub(r"[\s\*]", "", seq).upper()
    m = re.search("[^A-Z]", nseq)
    if m:
        raise RuntimeError("Normalized sequence contains non-alphabetic characters")
    return nseq
```

So the pipeline is: **strip `\s` and `*`** → **uppercase** → **reject
non-A-Z** → **SHA-512** → **truncate 24 bytes** → **base64url**.

### gtars RefgetStore

gtars is a Rust library with Python bindings. The `digest_sequence`
function (from `gtars/refget/__init__.pyi`, line 1663) documents:

> "The input data is automatically uppercased to ensure consistent digest
> computation (matching FASTA processing behavior)."

No mention of stripping `*` or any other characters. The Rust
implementation computes SHA-512 directly over the uppercased FASTA bytes,
including any `*` characters.

When `add_sequence_collection_from_fasta` ingests a FASTA file, each
sequence record's `sha512t24u` field is computed from the raw
(uppercased) FASTA content.

### Result

For a protein like ENSP00000340011.8 (297 aa, `*` at positions
239, 269, 280):

```
seqrepo seq_id:    fv4uQW4AfrwGj3dx9D4xfpccU7Qnqtia  ← hash of 294 chars (3 * stripped)
gtars sha512t24u:  8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0  ← hash of 297 chars (raw)
```

For any DNA sequence or any protein without `*`, the two algorithms
produce identical output.

## 4. What the RefGet specification says

The RefGet spec (v1.0.1 / v2.0) defines sequence normalization for
checksum computation:

> "When calculating the checksum for a sequence, all non-base symbols
> (\n, spaces, etc) must be removed and then uppercase the rest. The
> allowed alphabet for checksum calculation is uppercase ASCII
> (`0x41`-`0x5A` or `A-Z`)."

The `*` character (ASCII `0x2A`) is not in the range `0x41`–`0x5A`.
By the allowed-alphabet clause, `*` is a non-base symbol that should
be removed before hashing.

The VRS specification's type prefix table
([computed_identifiers.rst](https://github.com/ga4gh/vrs/blob/main/docs/source/conventions/computed_identifiers.rst))
delegates `SQ.*` identifier computation to RefGet:

```
SQ, Sequence (RefGet)
```

By the letter of the RefGet spec, **bioutils/seqrepo is conformant**
(strips `*` before hashing) and **gtars is non-conformant** (hashes
`*` with the sequence).

### Ambiguity

The spec says "non-base symbols (\n, spaces, **etc**)" — the "etc" is
vague, and `*` is a legitimate IUPAC amino acid code (Ter, adopted in
HGVS nomenclature v2.0). The spec was designed for nucleotide sequences
and says of proteins only: "other sequences could be provided via the
same mechanism e.g. cDNA, CDS, mRNA or proteins." There is no
protein-specific normalization guidance.

A reasonable case can be made either way:
- **Strip**: `*` is not in A-Z and the spec says A-Z only
- **Keep**: `*` is a legitimate IUPAC code encoding real positional
  information, and the "non-base symbols" examples suggest formatting
  characters (whitespace, line breaks), not sequence content

### Collision risk from stripping

Stripping `*` creates deliberate collisions: `"AB*CD"` and `"ABCD"`
produce the same digest. In practice this risk is negligible — these
proteins don't have `*`-free counterparts in any real database. None of
the 115 affected ENSPs shares a seqrepo `seq_id` with any NP_ or XP_
accession (verified: 0/82 unique seq_ids overlap with NCBI namespace).

## 5. Why this was invisible before Ensembl

NCBI RefSeq protein FASTAs (`human.*.protein.faa.gz`) contain zero `*`
characters across all 197,939 sequences in all 16 shards:

```
RefSeq protein shards (16 files):  197,939 sequences, 0 lines with *
Ensembl pep FASTA (1 file):        123,845 sequences, 158 with * (115 in seqrepo)
```

NCBI's convention: if translation hits a stop codon, no protein FASTA is
generated. Ensembl's convention: annotate the CDS anyway and emit `*` in
the protein FASTA. Since the `*`-stripping code path in bioutils only
changes the digest when `*` is actually present, the entire RefSeq
protein namespace (and all DNA sequences) produce identical digests in
both systems.

## 6. Ensembl source files for reproduction

The Ensembl protein FASTA that triggers this discrepancy is:

```
Homo_sapiens.GRCh38.pep.all.fa.gz
```

Downloaded from: `ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz`

This is the only file needed to reproduce the issue in a minimal test.
It is configured as a `[[seqset]]` entry in `misc/refgetstore/assemblies.toml`.

To reproduce, load this single FASTA into both a seqrepo instance and a
gtars RefgetStore, then compare the digests for the accessions listed
below.

## 7. Full list of affected accessions

116 total: 115 `star_normalization` (identical bytes, different digests)
+ 1 `sequence_divergence` (different bytes across Ensembl releases).

| Accession | Length | Stars | Gene | Biotype | Cause | SeqRepo seq_id | gtars sha512t24u |
|---|---|---|---|---|---|---|---|
| ENSP00000340011.8 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | fv4uQW4AfrwGj3dx9D4xfpccU7Qnqtia | 8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0 |
| ENSP00000375609.1 | 225 | 1 | KIR2DS4 | protein_coding | star_normalization | 0wFdmvBugnimf4tOFjCVdckpw3P6GmOx | Qg7FS7qyO-pfHp7HYd4wvr26DEx1oov4 |
| ENSP00000385765.6 | 230 | 1 | GSTT2 | protein_coding | star_normalization | lV7Y-VqQL6dNPJ2ZKWTqw46p_da8Eodm | YL07fWvAjY43qO_0ZnVw5BxCvo_M_0x- |
| ENSP00000402491.1 | 205 | 1 | AKR7L | protein_coding | star_normalization | e6sQaY-e_N_HsVyuY3ymiFysYBMBFk9W | dBxeQN4GInX8eFSsicpYD8N5Nn_0My2g |
| ENSP00000411913.1 | 331 | 1 | AKR7L | protein_coding | star_normalization | _RbUr7Oq9qE22xU9fnMC2uLA8RIp2Hri | TQ93TCrFm5e-mb-tkSHz3d_DYaAlo4i3 |
| ENSP00000421815.1 | 226 | 1 | CASP12 | protein_coding | star_normalization | tM9WHP2MFnlP7TaVWy-zF3pASsoLzDdU | PuFuh2fJaeBhwf_CIXb9c4n-5HHCkIA8 |
| ENSP00000422220.1 | 162 | 1 | GBA3 | protein_coding | star_normalization | o6IopMWmnKAE7FUNdEFy2TnsyBDa-swp | 2ludFB9wsyqzj2rD-ba2o2sFeXQIa277 |
| ENSP00000423754.1 | 298 | 1 | GBA3 | protein_coding | star_normalization | K56fM_BMAirMsSULhhLEc2yBkgvEE_aZ | kw3yUvbNK6RIftnKndZqroj4_OZxeNl9 |
| ENSP00000423899.1 | 224 | 1 | CASP12 | protein_coding | star_normalization | RDSLXXiECCksD3fd4HEjkDRDUojB3zm1 | D9gmmANCDpy0XX2JYSF7ZeHkXVaaDAm6 |
| ENSP00000423970.1 | 295 | 1 | CASP12 | protein_coding | star_normalization | fOlJtYkBTNNq6H0kRnKsGbSqtO5PcFEr | 96AGBZz3OkRXlz54m0a5vf04XVfn8MG0 |
| ENSP00000424038.1 | 341 | 1 | CASP12 | protein_coding | star_normalization | lOkFBmBx4LkR4GtP4ljxF4ssKcMpU9S2 | sNGrVi9jJokpvjhGaAkw9-5mCjrkxR5s |
| ENSP00000425652.1 | 277 | 1 | CASP12 | protein_coding | star_normalization | 0UyEhk9QGqb2gKNUMpHZ93nXVAXK9yTO | JQU8ZqYvYRuVO19ZDaHQhHM89PyidybZ |
| ENSP00000426427.1 | 140 | 1 | CASP12 | protein_coding | star_normalization | 03nXBcwbbJQPYE_5LurPSoCpnsikUdEn | pO53MQMJFDxQT6bOb7ey9est-55eH5NH |
| ENSP00000426566.1 | 257 | 1 | CASP12 | protein_coding | star_normalization | _nkNH8COTskYyNkCoTzIRNHzENVylj1y | egHuRL5C7VrWvOx3FvqbEayoBNIeC6ZM |
| ENSP00000426627.1 | 323 | 1 | FCGR2C | protein_coding | star_normalization | OpwKkqOENTdevZ_CrLEOVumim2QTctaK | 2KlIfbXxNKZoNImdwrmKBdRYXPkKL_2G |
| ENSP00000427437.1 | 142 | 1 | CASP12 | protein_coding | star_normalization | wWUGUKB2gnrfgmMm0qg34Iy9mUEhx4ju | 7KGhHMvva3lM9749_TxBRCl5Oz8YcdKo |
| ENSP00000427458.1 | 469 | 1 | GBA3 | protein_coding | star_normalization | rFjG0Bt4bfEH09x5PaR9jZGd03wmgE6O | -x9lIom6XRVzSbgQDyNb0hrandBbttoV |
| ENSP00000427969.1 | 10 | 2 | IGHD2-15 | IG_D_gene | star_normalization | y9mM4wkXZleVNSyYRPm83MhQd47vpV0A | Ah8IodVoXqfHk4DDt8jsFGSic0fRZYU0 |
| ENSP00000428366.1 | 12 | 1 | IGHD3-16 | IG_D_gene | star_normalization | SHyk_qyvZUk1fYzNdSTs6Drv8jS87r_C | CN6dOnzsnNZiifsD2mNJ941p8KVRVPJu |
| ENSP00000428393.1 | 5 | 2 | IGHD4-4 | IG_D_gene | star_normalization | qRfoOZ40q93mvxVq3MDnu2rMZt2mB_ob | yAZ-3aLns9fo2SxCNtvZJGjw_9-p67nz |
| ENSP00000428616.1 | 10 | 1 | IGHD2-8 | IG_D_gene | star_normalization | KRLK6ZBMvRpTakgKIIWxTzvHu6nZtEKQ | 5ZDMzPF5f4u5xiuBDA2ByfCYoG4534KM |
| ENSP00000429324.1 | 9 | 1 | IGHD2-21 | IG_D_gene | star_normalization | m19OoHXspP1j74FYLLR0Fon_Tca6IwAP | fn4oRtmkRm5ikN4QPXXg0RdDOzaLpvUH |
| ENSP00000429952.1 | 10 | 3 | IGHD3-22 | IG_D_gene | star_normalization | 5Ek4xsz2ZS7GfXgv_ZBjOH5FkAP3XMX6 | JCCHtapvgLCCSu0eIJCK34DXB_mOSDrc |
| ENSP00000430034.1 | 5 | 2 | IGHD4-11 | IG_D_gene | star_normalization | qRfoOZ40q93mvxVq3MDnu2rMZt2mB_ob | yAZ-3aLns9fo2SxCNtvZJGjw_9-p67nz |
| ENSP00000430248.1 | 6 | 2 | IGHD4-23 | IG_D_gene | star_normalization | AAWZIb3sgZG4DGzYKPTEdsDRPppA20nl | ZyuIuWLGyac0hNqMfiw0XIJIjTx1UF2Y |
| ENSP00000430788.1 | 10 | 2 | IGHD2-2 | IG_D_gene | star_normalization | 4nJC7z8ekuHU6kAl4yMiV6qZlrpomr_6 | vQGG6MQqX_mw1KxS82QU0Kuk5KtUItWA |
| ENSP00000431031.1 | 1318 | 1 | MROH5 | protein_coding | star_normalization | l4SXGqAbFUh7r1jc-ImcNEdp067Fa0MY | zkX9kDvYwL768RXctEHW7-3lneBHc71d |
| ENSP00000431089.1 | 5 | 2 | IGHD4-17 | IG_D_gene | star_normalization | BM1Hbmp96w4OSkWZ5u_ROCq4tWMf_V69 | hzlOSDnnFBtjxI-qUTYSXVvP8K76MfHk |
| ENSP00000434417.1 | 2459 | 17 | PKD1L2 | protein_coding | star_normalization | zKGzcqTwIi1cpsh_25jTR-g7cycsIEpf | Ih9FQHBcIpxLvdvSmQwf1B_uW07yCbjE |
| ENSP00000434644.1 | 1774 | 17 | PKD1L2 | protein_coding | star_normalization | DZRDVc8gAxzoOkjyLBLK3IcVRKKK46Pw | v8czIuEX2XTC0M9Dw7Ax3fzb4H8qdgOL |
| ENSP00000438013.1 | 314 | 4 | OR4Q2 | protein_coding | star_normalization | nDmCYyT1EKuf9Z8mrnEgCxgF2Hjn0T_P | _nPal56b66hTKeKlOajKFSKjKY8DcuU- |
| ENSP00000439604.1 | 516 | 11 | CYP2D7 | protein_coding | star_normalization | oaVc-acao8fOMIyVp_H98dcz0BSDIUWH | aYF_Y_ycp3u3t9Gqeu96GLzc304fQ0zp |
| ENSP00000445124.1 | 497 | 11 | CYP2D7 | protein_coding | star_normalization | 6VmdTL7nqa7JQpolg-wuBLGwIKBnrP9F | MIwe_aAViS3zDJuv3mQTnIDVQMSLST8L |
| ENSP00000451021.2 | 314 | 1 | OR11H7 | protein_coding | star_normalization | 38LKEEE-120Kwsih-vNlwnx8tb2A4ZSx | _rwlFqSN0jljdGdyuRQXsthR-SSMBJ2G |
| ENSP00000452150.1 | 21 | 1 | TRAJ58 | TR_J_gene | star_normalization | ynPhtGE3u23nkLJNMxVv51pviMYUrBBX | Ede4Q4UXqG9Dzi4MrnjeH9jp9k0ggu7T |
| ENSP00000461363.1 | 517 | 11 | CYP2D7 | protein_coding | star_normalization | O0rR6w_qtUdlVU360MFqtPEbi7ZFraqF | KInA0nFf_sQ_yQXjq92_k4rkC72kfgDG |
| ENSP00000461768.1 | 498 | 11 | CYP2D7 | protein_coding | star_normalization | SBAWxYVz_UK9DZodtRFfbK1-2l5-uvL2 | 6UVE7AZ23iHnnV0dDaNAFu_wGKHPYYoI |
| ENSP00000463502.4 | 469 | 1 | PNLIPRP2 | protein_coding | star_normalization | syXOnAwOvxwpfLO67Z7c84KArOxWUGX- | zWPI0igbczAoDoC3RNroPzMKmMmNOe6I |
| ENSP00000474017.2 | 10 | 2 | IGHD2OR15-2B | IG_D_gene | star_normalization | PkCr3cRIsERXzmDHBNwJKKU2AhW1pcjS | a8CMClpZG-W05aaW-ZTFsdZ8WUQnX6t4 |
| ENSP00000474065.2 | 10 | 2 | IGHD2OR15-2A | IG_D_gene | star_normalization | PkCr3cRIsERXzmDHBNwJKKU2AhW1pcjS | a8CMClpZG-W05aaW-ZTFsdZ8WUQnX6t4 |
| ENSP00000474133.2 | 10 | 1 | IGHD3OR15-3A | IG_D_gene | star_normalization | U7neqpeNT4iNk-l3yhxlfzeh7imlPahR | Qflup9oKLA9_ebzhyw1VwhdIcInZGnb_ |
| ENSP00000474573.2 | 10 | 1 | IGHD3OR15-3B | IG_D_gene | star_normalization | U7neqpeNT4iNk-l3yhxlfzeh7imlPahR | Qflup9oKLA9_ebzhyw1VwhdIcInZGnb_ |
| ENSP00000474693.2 | 6 | 2 | IGHD4OR15-4A | IG_D_gene | star_normalization | FfJUWv8y8zqeEkxc8WC0CBxAVBpHKby_ | SMi7IfQFm46hxfH46aUlNOJd_mTj6CYC |
| ENSP00000475053.2 | 6 | 2 | IGHD4OR15-4B | IG_D_gene | star_normalization | FfJUWv8y8zqeEkxc8WC0CBxAVBpHKby_ | SMi7IfQFm46hxfH46aUlNOJd_mTj6CYC |
| ENSP00000475160.1 | 308 | 9 | OR5AC1 | protein_coding | star_normalization | FOZR2eNkiqv0vQKW28sYDoFUk-pHLf7q | GO28PhtqjbkgqI_ZXen46QVdAUtQsTy0 |
| ENSP00000475351.1 | 315 | 8 | OR12D1 | protein_coding | star_normalization | NaBv12DVCEzFhBS1k7fnmH1RvEOocOQZ | JcuCvOWNKvugDHiBlc5QNxiD9c1xGoHK |
| ENSP00000475611.1 | 311 | 6 | OR10J4 | protein_coding | star_normalization | CafOAbOr7B4OwynxW2mGBs93mmWzlLpj | so4CHY-xDA2Kza_d9Sq8Cz9E96yZMrL8 |
| ENSP00000476186.1 | 310 | 2 | OR5H8 | protein_coding | star_normalization | GQpNKfcojbpwHYXfeaYzPUelZhAsVouD | Mg_qtBntkSHm5xPGBWSLyzHBXF1Ua1Oi |
| ENSP00000476380.1 | 326 | 1 | OR10AC1 | protein_coding | star_normalization | 546P9zkeSm9jcCq9v5rNou3-U-9wOXa0 | 9gXvvtU-EONZkDaTQZTRkqqAeX9YbCjD |
| ENSP00000476467.1 | 315 | 1 | OR4A8 | protein_coding | star_normalization | HV7AX9iHZz2qgacGsUkeLl9JY1wQSy7b | 91F6s9En0f8Tp5H18TUHjtFIhnFA79lA |
| ENSP00000476576.1 | 308 | 3 | OR52E1 | protein_coding | star_normalization | MHCBdS-C2GvAKyaWy5g7YpM5Tl6S7r2P | exEjrMOOllerXl47soSNdLchtvUnq-n4 |
| ENSP00000477333.1 | 319 | 12 | OR13C7 | protein_coding | star_normalization | 81G0gd-_pNefpEX9qCl6Xo3SfpL_GtM_ | jBGqsI0x-vXuNDpnaGRFodo1JLeIrtqk |
| ENSP00000477337.2 | 314 | 16 | OR5G3 | protein_coding | star_normalization | Hgx6zf3aYxje2c12BUxwnW3Wfwcz72-E | sfkGeZ9IhJISVjv1K1tsfeX8obIqTnf0 |
| ENSP00000477393.1 | 328 | 10 | OR5AL1 | protein_coding | star_normalization | -a7aCjSuO0y-Noz6SjHENVXwfOXx2lgl | C8md75ie9b7ZYfCXTAlyTnIUF4YkeQkn |
| ENSP00000478383.1 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | se_guVGDbGUwQ_VwGptcBPiTSUEMOTDA | IlMIcqDnTyx9LM1oe18o-6-Q6QWvov_6 |
| ENSP00000478656.2 | 517 | 11 | CYP2D7 | protein_coding | star_normalization | 8Jwk52z5pfAcprdAXS7iQN1kTNxxio40 | EbKxx4ZLK27lrDbLCoPaQMFIuY1_aMQs |
| ENSP00000478998.1 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | Ot3JN_E9VBadQ46vgolO3JMynJDqSGHa | f01vgyglFX1ju7emFuoB3OQn6mBWjCR6 |
| ENSP00000479129.1 | 225 | 1 | KIR2DS4 | protein_coding | star_normalization | PVMp-pPMZHSdjxIHbf9BPZFTVliQPwt8 | TiIdLvCGnMqH2aeTyjlo9IlWOwojBX4g |
| ENSP00000479356.1 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | se_guVGDbGUwQ_VwGptcBPiTSUEMOTDA | IlMIcqDnTyx9LM1oe18o-6-Q6QWvov_6 |
| ENSP00000479369.1 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | fv4uQW4AfrwGj3dx9D4xfpccU7Qnqtia | 8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0 |
| ENSP00000479404.1 | 374 | 4 | KIR2DS4 | protein_coding | star_normalization | t3g70NPf5GGSyB223lAwFe6zw1ChtmEy | F_y21mXr_NuyB93ooYdRBy0kd-PXP6tz |
| ENSP00000481139.1 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | 6y4nN8mld_LcvD-EbegjezG4J0ipyLq1 | GIq_KqgXV64sQO3gTLidkdIOtqkmWH7A |
| ENSP00000481279.1 | 374 | 4 | KIR2DS4 | protein_coding | star_normalization | XwjGxQ5MlPAw2PL52E6hJdBXJE3Z-Q32 | cUaEhxJnZ5QAuY3JG0Bn6bcW9MwKq73n |
| ENSP00000481428.1 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | fv4uQW4AfrwGj3dx9D4xfpccU7Qnqtia | 8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0 |
| ENSP00000482055.2 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | se_guVGDbGUwQ_VwGptcBPiTSUEMOTDA | IlMIcqDnTyx9LM1oe18o-6-Q6QWvov_6 |
| ENSP00000482745.1 | 341 | 1 | CASP12 | protein_coding | star_normalization | lOkFBmBx4LkR4GtP4ljxF4ssKcMpU9S2 | sNGrVi9jJokpvjhGaAkw9-5mCjrkxR5s |
| ENSP00000483571.1 | 225 | 1 | KIR2DS4 | protein_coding | star_normalization | lv5JmDOztbvJgKGluw0Fh6jFTb0LyBFT | mNBJqEx1npRp6Qy7ofSZAYW4iX81-b3m |
| ENSP00000483620.1 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | fv4uQW4AfrwGj3dx9D4xfpccU7Qnqtia | 8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0 |
| ENSP00000485881.1 | 311 | 15 | OR5G3 | protein_coding | star_normalization | 5CcawTtWZa2Bc5otcdnSdzOQV1aMT_cq | N3f7prhPB-ztQHfhGD48lojOdU42AWwG |
| ENSP00000485975.1 | 315 | 1 | OR8J2 | protein_coding | star_normalization | AlKSb7ywKLef2uY_oqJ97OqTeZ1Cdb1D | SzZa7T9y81S_V8iQp9RnpbnKshXX9IHp |
| ENSP00000486295.1 | 4040 | 11 | MUC6 | protein_coding | star_normalization | iyPvjohGDgVYeFEVYQ_n6DFS5sOkwUTI | VJMwcSEKzmAksrvE276Okup2ueBujczb |
| ENSP00000487044.1 | 498 | 11 | CYP2D7 | protein_coding | star_normalization | HeRZbi709aJqVjELYnnZZJrRUW17WjWH | mdGiXxUvar3WuR1h-zsoryFcpLHQgi0H |
| ENSP00000487059.1 | 4714 | 23 | MUC6 | protein_coding | star_normalization | NpQlnZ0jUAHvbmGF37prMDgMPv3THHYN | 2Vgy0wYjqGFzUpqGWPDJipU0XNdrKU8a |
| ENSP00000487596.1 | 493 | 1 | MROH5 | protein_coding | star_normalization | _IPtIdLGP1_WEK4al_dKDId2NigEHXZn | bZJHWz7sSJQUC00fzSqS0wgZMcQ4X1i4 |
| ENSP00000487599.1 | 9 | 1 | IGHD2-21 | IG_D_gene | star_normalization | m19OoHXspP1j74FYLLR0Fon_Tca6IwAP | fn4oRtmkRm5ikN4QPXXg0RdDOzaLpvUH |
| ENSP00000487604.1 | 10 | 2 | IGHD2-2 | IG_D_gene | star_normalization | 4nJC7z8ekuHU6kAl4yMiV6qZlrpomr_6 | vQGG6MQqX_mw1KxS82QU0Kuk5KtUItWA |
| ENSP00000487787.1 | 12 | 1 | IGHD3-16 | IG_D_gene | star_normalization | SHyk_qyvZUk1fYzNdSTs6Drv8jS87r_C | CN6dOnzsnNZiifsD2mNJ941p8KVRVPJu |
| ENSP00000487993.1 | 10 | 2 | IGHD2-15 | IG_D_gene | star_normalization | y9mM4wkXZleVNSyYRPm83MhQd47vpV0A | Ah8IodVoXqfHk4DDt8jsFGSic0fRZYU0 |
| ENSP00000488083.1 | 10 | 1 | IGHD2-8 | IG_D_gene | star_normalization | KRLK6ZBMvRpTakgKIIWxTzvHu6nZtEKQ | 5ZDMzPF5f4u5xiuBDA2ByfCYoG4534KM |
| ENSP00000488168.1 | 6 | 2 | IGHD4-23 | IG_D_gene | star_normalization | AAWZIb3sgZG4DGzYKPTEdsDRPppA20nl | ZyuIuWLGyac0hNqMfiw0XIJIjTx1UF2Y |
| ENSP00000488261.1 | 5 | 2 | IGHD4-17 | IG_D_gene | star_normalization | BM1Hbmp96w4OSkWZ5u_ROCq4tWMf_V69 | hzlOSDnnFBtjxI-qUTYSXVvP8K76MfHk |
| ENSP00000488314.1 | 153 | 1 | UBE2NL | protein_coding | star_normalization | Zbnpiff-E255gCpMrQGZejm0OauEumnc | 4yrFdP8maHkqUKDKBSVQ-gIrkehiwHxr |
| ENSP00000488711.1 | 10 | 3 | IGHD3-22 | IG_D_gene | star_normalization | 5Ek4xsz2ZS7GfXgv_ZBjOH5FkAP3XMX6 | JCCHtapvgLCCSu0eIJCK34DXB_mOSDrc |
| ENSP00000488735.1 | 5 | 2 | IGHD4-11 | IG_D_gene | star_normalization | qRfoOZ40q93mvxVq3MDnu2rMZt2mB_ob | yAZ-3aLns9fo2SxCNtvZJGjw_9-p67nz |
| ENSP00000488889.1 | 5 | 2 | IGHD4-4 | IG_D_gene | star_normalization | qRfoOZ40q93mvxVq3MDnu2rMZt2mB_ob | yAZ-3aLns9fo2SxCNtvZJGjw_9-p67nz |
| ENSP00000488911.1 | 497 | 12 | CYP2D7 | protein_coding | star_normalization | h9XYkJfTPZZjGI3uq_VgHAxQfGOnusMx | Lim-8--cZivPCz0LQqISF0tUGv2hrbrX |
| ENSP00000488993.1 | 244 | 1 | GSTT2 | protein_coding | star_normalization | Pu86r-b0I3U_dqkK9zQ2OPGSshhhhHWG | DnlfpzUcYbk15_dKw1jqPVCwqxRENcFT |
| ENSP00000489014.1 | 516 | 12 | CYP2D7 | protein_coding | star_normalization | -fkYMovNew_mTufEcBSEn80kG2rr4OB_ | ZurDzxk4uhTHtf8xZAPlkm2VY4kZwv5E |
| ENSP00000489240.1 | 108 | 1 | IFNL4 | protein_coding | star_normalization | aPj392cUfAWEeDcwOZLZafkPAwlUuscd | IoKaC3TiSvxPlN1u_b9QaSEGIr9uavgu |
| ENSP00000489559.1 | 132 | 1 | IFNL4 | protein_coding | star_normalization | 7Dqau58VfarSEwoAyKydRrGg8qHMU6sN | jdWrXmbGQ39RS_FxE9AxivWfxAjk0GgM |
| ENSP00000490976.1 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | Ot3JN_E9VBadQ46vgolO3JMynJDqSGHa | f01vgyglFX1ju7emFuoB3OQn6mBWjCR6 |
| ENSP00000491125.1 | 357 | 1 | VN1R5 | protein_coding | star_normalization | tRRBvR9lKL2Pqn0LDyPFWX9QzwN6R9cI | yJ1li88C4ZhzDqpmIMuJhDlayDIu2i15 |
| ENSP00000491133.1 | 306 | 1 | OR4C45 | protein_coding | star_normalization | 84ZbQcviHt30Ks_HfJnC25WT1H3Ll_aG | UPwfn83YcGH13BQeW32xlO0MSD6SS24V |
| ENSP00000491451.1 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | se_guVGDbGUwQ_VwGptcBPiTSUEMOTDA | IlMIcqDnTyx9LM1oe18o-6-Q6QWvov_6 |
| ENSP00000491543.1 | 225 | 1 | KIR2DS4 | protein_coding | star_normalization | Ipbqtow-myKtzhLMr1FiU01rBRs7Iwv1 | w6SG0sotavd93SmwWs3HxDb31icvb43V |
| ENSP00000491974.1 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | fv4uQW4AfrwGj3dx9D4xfpccU7Qnqtia | 8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0 |
| ENSP00000491994.1 | 225 | 1 | KIR2DS4 | protein_coding | star_normalization | 0wFdmvBugnimf4tOFjCVdckpw3P6GmOx | Qg7FS7qyO-pfHp7HYd4wvr26DEx1oov4 |
| ENSP00000492028.2 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | fv4uQW4AfrwGj3dx9D4xfpccU7Qnqtia | 8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0 |
| ENSP00000492090.1 | 225 | 1 | KIR2DS4 | protein_coding | star_normalization | Ipbqtow-myKtzhLMr1FiU01rBRs7Iwv1 | w6SG0sotavd93SmwWs3HxDb31icvb43V |
| ENSP00000492193.1 | 688 | 1 | POTEA | protein_coding | star_normalization | 7Pr4ISih0_g-We99919kiPT8Ft1voOjU | QKLfHnF-uMJbWpYlzxPWivnK9FISFaHv |
| ENSP00000492245.1 | 225 | 1 | KIR2DS4 | protein_coding | star_normalization | Ipbqtow-myKtzhLMr1FiU01rBRs7Iwv1 | w6SG0sotavd93SmwWs3HxDb31icvb43V |
| ENSP00000492265.1 | 642 | 1 | POTEA | protein_coding | star_normalization | NdYD2MxngnXGXKck56SC5Y2pwnSRxQeN | gZoUTjOaMrKyFO7XcjzV7e0F6hHw0Y2I |
| ENSP00000492337.2 | 225 | 1 | KIR2DS4 | protein_coding | star_normalization | 0wFdmvBugnimf4tOFjCVdckpw3P6GmOx | Qg7FS7qyO-pfHp7HYd4wvr26DEx1oov4 |
| ENSP00000492458.1 | 225 | 1 | KIR2DS4 | protein_coding | star_normalization | PVMp-pPMZHSdjxIHbf9BPZFTVliQPwt8 | TiIdLvCGnMqH2aeTyjlo9IlWOwojBX4g |
| ENSP00000492627.1 | 297 | 3 | KIR2DS4 | protein_coding | star_normalization | se_guVGDbGUwQ_VwGptcBPiTSUEMOTDA | IlMIcqDnTyx9LM1oe18o-6-Q6QWvov_6 |
| ENSP00000493090.1 | 331 | 2 | OR1P1 | protein_coding | star_normalization | OwOpiUeTZxOiRyJpUVgFPEknheNaytQk | 7gOPz5rjQBMvXPqGeF8kBFb1MRWsv-XP |
| ENSP00000493099.2 | 105 | 1 | SCYGR10 | protein_coding | star_normalization | lcX6lIyYgyvgNIcxL8Nf_9c2te0uWR7J | bX-m0VrO77yaMLzVcz__B14CYNhMhNS8 |
| ENSP00000493221.1 | 317 | 6 | OR2T7 | protein_coding | star_normalization | YhrPozOQOjywrkWpQWPgu49GdiZqFZk8 | vax-AZkewD5SpDAT5szm8meLHgf2beXC |
| ENSP00000493243.1 | 317 | 6 | OR2T7 | protein_coding | star_normalization | YhrPozOQOjywrkWpQWPgu49GdiZqFZk8 | vax-AZkewD5SpDAT5szm8meLHgf2beXC |
| ENSP00000493291.1 | 315 | 7 | OR4K3 | protein_coding | star_normalization | _CbilrTwIJnxQ8UDDZ3ieqBW5Jz6IlVd | 8GVt4hPQBAH0n47yoEBplEk176LMj6fZ |
| ENSP00000493317.1 | 315 | 1 | OR8J2 | protein_coding | star_normalization | AlKSb7ywKLef2uY_oqJ97OqTeZ1Cdb1D | SzZa7T9y81S_V8iQp9RnpbnKshXX9IHp |
| ENSP00000493452.1 | 308 | 1 | OR9H1P | protein_coding | star_normalization | mICKQK_3Sn5lQNys1NRoohE3GrKeB5eh | HO0I-zZFjRZAJMqZ7io8-u7lxa3y9P6_ |
| ENSP00000495205.1 | 954 | 1 | PTCHD3 | protein_coding | star_normalization | Nmuo8TKv6WCrdaDDfamxvpgrA26KLpLS | qR9P0vEzxBhBS34XQHLkcbx9GCG08qoa |
| ENSP00000496497.1 | 326 | 1 | OR10AC1 | protein_coding | star_normalization | 546P9zkeSm9jcCq9v5rNou3-U-9wOXa0 | 9gXvvtU-EONZkDaTQZTRkqqAeX9YbCjD |
| ENSP00000498832.1 | 232 | 2 | FAM246C | protein_coding | star_normalization | PZVNRQLqkqRCjj4klpDJmfxJqLHzCSyO | Gr7F6_lh_SAfec-UJ4tgY0vLz7oag0LV |
| ENST00000668831.1 | 1120 | 0 | | | sequence_divergence | NllKO15Wb60jTlfVHY5YOFfN3WFZyn5M | 2EdYgXlBf5xGHlhCifMYtWt9QZPGOyxA |

## 8. References

### Specifications

- [RefGet Sequences v2.0](https://ga4gh.github.io/refget/sequences/) — "allowed alphabet for checksum calculation is uppercase ASCII (0x41-0x5A or A-Z)"
- [RefGet v1.0.1 (hts-specs)](https://samtools.github.io/hts-specs/refget.html) — "all non-base symbols (\n, spaces, etc) must be removed"
- [RefGet spec source (GitHub)](https://github.com/samtools/hts-specs/blob/master/refget.md)
- [VRS Computed Identifiers](https://github.com/ga4gh/vrs/blob/main/docs/source/conventions/computed_identifiers.rst) — SQ type prefix delegates to RefGet
- [HGVS amino acid nomenclature](https://www.hgvs.org/mutnomen/references.html) — `*` / Ter as IUPAC code for stop

### Source code

- `bioutils.sequences.normalize_sequence` — `.venv/.../bioutils/sequences.py` — `re.sub(r"[\s\*]", "", seq).upper()`
- `bioutils.digests.seq_seqhash` — `.venv/.../bioutils/digests.py` — calls `normalize_sequence` then `vmc_digest`
- `biocommons.seqrepo.seqrepo.SeqRepo` — `.venv/.../biocommons/seqrepo/seqrepo.py` — line 11: `from bioutils.digests import seq_seqhash as sha512t24u`; line 148: `seqhash = sha512t24u(seq)` (actually `seq_seqhash`)
- `gtars.refget.digest_sequence` — Rust implementation, Python stub at `.venv/.../gtars/refget/__init__.pyi:1663` — "input data is automatically uppercased" (no `*` stripping)
- `gtars.refget.sha512t24u_digest` — `.venv/.../gtars/refget/__init__.pyi:1544` — raw hash, no normalization

### Diagnostic artifacts

- `misc/refgetstore/diagnose_ensembl_mismatches.py` — 4-phase diagnostic script
- `misc/refgetstore/diagnosis_report.json` — full JSON report
- `misc/refgetstore/ensembl_known_divergent.txt` — TSV of affected accessions
