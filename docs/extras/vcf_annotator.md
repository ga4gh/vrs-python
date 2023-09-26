# VCF Annotator

The [VCF Annotator tool](../../src/ga4gh/vrs/extras/vcf_annotation.py) provides utility for annotating VCF's with VRS Allele IDs.

## How to use

*Note:\
The examples run from the root of the vrs-python directory and assumes that `input.vcf.gz` lives in the current directory*

To see the help page:

```commandline
python3 -m src.ga4gh.vrs.extras.vcf_annotation --help
```

### Use local SeqRepo Data Proxy with default root directory

The tool uses a SeqRepo data proxy. By default, the local instance at `/usr/local/share/seqrepo/latest` is used.

Example of how to run:

```commandline
python3 -m src.ga4gh.vrs.extras.vcf_annotation --vcf_in input.vcf.gz --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl
```

`--vcf_in` specifies the path of the input VCF file to annotate. `--vcf_out` specifies the path of the output annotated VCF file. The `--vrs_pickle_out` specifies the path of the output pickle file containing VRS data.

### Use local SeqRepo Data Proxy with different

You can change the root directory of SeqRepo by using `seqrepo_root_dir`.

To use the local SeqRepo data proxy with SeqRepo root directory at `vrs-python/seqrepo/latest`:

```commandline
python3 -m src.ga4gh.vrs.extras.vcf_annotation --vcf_in input.vcf.gz --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl --seqrepo_root_dir vrs-python/seqrepo/latest
```

### Use the REST SeqRepo Data Proxy with default base url

You can change the data proxy type by using: `--seqrepo_dp_type` (options are `local` or `rest`).

To use the REST SeqRepo data proxy at default url: `http://localhost:5000/seqrepo`:

```commandline
python3 -m src.ga4gh.vrs.extras.vcf_annotation --vcf_in input.vcf.gz --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl --seqrepo_dp_type rest
```

### Use the REST SeqRepo Data Proxy with different base url
You can change the SeqRepo REST base url by using: `--seqrepo_base_url`.

To use the REST SeqRepo data proxy, at custom url: `http://custom.url:5000/seqrepo`:
```commandline
python3 -m src.ga4gh.vrs.extras.vcf_annotation --vcf_in input.vcf.gz --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl --seqrepo_dp_type rest --seqrepo_base_url http://custom.url:5000/seqrepo
```
