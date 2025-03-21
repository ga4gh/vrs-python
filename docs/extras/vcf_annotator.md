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

Like other VRS-Python tools, the VCF annotator requires access to [sequence and identifier data services](https://vrs.ga4gh.org/en/stable/impl-guide/required_data.html#data-services), as implemented in libraries like [SeqRepo](https://github.com/biocommons/biocommons.seqrepo). By default, the CLI will attempt to connect to a [SeqRepo REST instance](https://github.com/biocommons/seqrepo-rest-service) at `http://localhost:5000/seqrepo`, but a URI can be passed with the `--dataproxy_uri` option or set with the `GA4GH_VRS_DATAPROXY_URI` environment variable (the former takes priority over the latter).

For example, to use a local set of SeqRepo data, you can use an absolute file path:

```commandline
vrs-annotate vcf --dataproxy_uri="seqrepo+file:///usr/local/share/seqrepo/2024-12-20/" --vcf_out=out.vcf.gz input.vcf.gz
```

Alternative, a relative file path:

```commandline
vrs-annotate vcf --dataproxy_uri="seqrepo+../seqrepo/2024-12-20/" --vcf_out=out.vcf.gz input.vcf.gz
```

Or an alternate REST path:

```commandline
vrs-annotate vcf --dataproxy_uri="seqrepo+http://mylabwebsite.org/seqrepo" --vcf_out=out.vcf.gz input.vcf.gz
```

### Other Options
`--vrs_attributes`
>Will include VRS_Start, VRS_End, VRS_State fields in the INFO field.

`--assembly` [TEXT]
>The assembly that the `vcf_in` data uses. [default: GRCh38]

`--skip_ref`
>Skip VRS computation for REF alleles.

`--require_validation`
>Require validation checks to pass in order to return a VRS object

`--help`
>Show the options available
