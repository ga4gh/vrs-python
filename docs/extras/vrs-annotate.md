# VRS annotate

Use `vrs-annotate` to annotate input files with VRS variation objects:

```shell
vrs-annotate --help
```

## VCF

Install `ga4gh.vrs[extras]` and configure [SeqRepo](../data-services.md).
Then annotate a VCF with VRS allele identifiers:

```shell
vrs-annotate vcf \
  --dataproxy-uri='seqrepo+http://localhost:5000/seqrepo' \
  --assembly=GRCh38 \
  --vcf-out=out.vcf.gz \
  input.vcf.gz
```

Set `--assembly` to the input VCF's reference assembly. Set the data-proxy URI with
`--dataproxy-uri` or `GA4GH_VRS_DATAPROXY_URI`. The command-line option overrides
the environment variable.

For available options and defaults:

```shell
vrs-annotate vcf --help
```

To annotate from Python, see the [VCF annotation API](../api/vcf-annotation.md).
