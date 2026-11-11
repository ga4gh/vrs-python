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

## Reading VRS Alleles from an annotated VCF

The [`VcfReader`](../../src/ga4gh/vrs/extras/annotator/vcf.py) class reconstructs the VRS Alleles from a VCF annotated with `--vrs-attributes` (or `vrs_attributes=True` in `VcfAnnotator.annotate()`). It iterates over the VCF's records and yields each pysam record together with a list of VRS Alleles:

```python
from pathlib import Path

from ga4gh.vrs.dataproxy import create_dataproxy
from ga4gh.vrs.extras.annotator.vcf import VcfReader

data_proxy = create_dataproxy("seqrepo+http://localhost:5000/seqrepo")
reader = VcfReader(data_proxy)
for record, vrs_alleles in reader.read(Path("out.vcf.gz"), assembly="GRCh38"):
    for vrs_allele in vrs_alleles:
        if vrs_allele is not None:
            print(record.chrom, record.pos, vrs_allele.id)
```

* The list of VRS Alleles has one value per REF and ALT allele (`record.alleles`), or per ALT allele (`record.alts`) if the VCF was annotated with `--skip-ref`. It holds `None` for alleles without a VRS ID, such as `*` ALTs, and for every allele of a record that failed annotation (see its `VRS_Error` INFO field).
* The data proxy is only used to look up the refget accession for each `CHROM`, and `assembly` must be the assembly used to annotate the VCF.
* The identifier of each reconstructed VRS Allele is checked against its `VRS_Allele_IDs` annotation. A `VcfReaderError` is raised if they differ, if the annotations of a record are incomplete, or if the VCF wasn't annotated with VRS attributes.
* `ReferenceLengthExpression` sequences longer than 50 residues aren't written to the annotated VCF, so these states are read without a `sequence`. This doesn't change their VRS identifiers.
