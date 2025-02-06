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

env var "GA4GH_VRS_DATAPROXY_URI" or CLI option

show all URI values

TODO

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
