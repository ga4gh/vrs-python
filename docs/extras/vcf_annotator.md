# VCF Annotator

The [VCF Annotator tool](../../src/ga4gh/vrs/extras/annotator/vcf.py) provides a Python class for annotating VCFs with VRS Allele IDs. A [command-line interface](../../src/ga4gh/vrs/extras/annotator/cli.py) is available for accessing these functions from a shell or shell script.

## How to use

*Note:\
The examples run from the root of the vrs-python directory and assumes that `input.vcf.gz` lives in the current directory*

To see the help page:

```commandline
vrs-annotate vcf --help
```

### Configuring the sequence data proxy

Like other VRS-Python tools, the VCF annotator requires access to [sequence and identifier data services](https://vrs.ga4gh.org/en/stable/impl-guide/required_data.html#data-services), as implemented in libraries like [SeqRepo](https://github.com/biocommons/biocommons.seqrepo). By default, the CLI will attempt to connect to a [SeqRepo REST instance](https://github.com/biocommons/seqrepo-rest-service) at `http://localhost:5000/seqrepo`, but a URI can be passed with the `--dataproxy-uri` option or set with the `GA4GH_VRS_DATAPROXY_URI` environment variable (the former takes priority over the latter).

For example, to use a local set of SeqRepo data, you can use an absolute file path:

```commandline
vrs-annotate vcf --dataproxy-uri="seqrepo+file:///usr/local/share/seqrepo/2024-12-20/" --vcf-out=out.vcf.gz input.vcf.gz
```

Alternative, a relative file path:

```commandline
vrs-annotate vcf --dataproxy-uri="seqrepo+../seqrepo/2024-12-20/" --vcf-out=out.vcf.gz input.vcf.gz
```

Or an alternate REST path:

```commandline
vrs-annotate vcf --dataproxy-uri="seqrepo+http://mylabwebsite.org/seqrepo" --vcf-out=out.vcf.gz input.vcf.gz
```

### Configuring optional RLE sequences

Set `GA4GH_VRS_RLE_SEQ_LIMIT` before starting Python or `vrs-annotate` to change
the maximum length of optional `ReferenceLengthExpression.sequence` values.
The default is `50` bases. Use a non-negative integer, `0` to omit the optional
sequence, or `none` (case-insensitive) for no limit. Invalid values raise a
`ValueError` when VRS-Python is imported.

The normalizer, translator, and VCF writer share this setting. Explicit
`rle_seq_limit` arguments to normalization or translation override the default;
the VCF writer still applies the configured output limit. Required
`LiteralSequenceExpression.sequence` values are unaffected.

For example:

```shell
GA4GH_VRS_RLE_SEQ_LIMIT=100 vrs-annotate vcf --vrs-attributes --vcf-out=out.vcf input.vcf
```

### Other Options
`--vrs-attributes`
>Will include VRS_Start, VRS_End, VRS_State fields in the INFO field.

`--assembly` [TEXT]
>The assembly that the `vcf-in` data uses. [default: GRCh38]

`--skip-ref`
>Skip VRS computation for REF alleles.

`--require-validation`
>Require validation checks to pass in order to return a VRS object

`--help`
>Show the options available
