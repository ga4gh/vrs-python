# Docker

This README includes information on how to use the docker containers.

## Running the VCF Annotator

Commands to run inside the Docker CLI. Ensure that the `vcf_out` and `vrs_file` are output to the `/data` directory so that data persists.

To use local SeqRepo (uses test data for example):

```commandline
python3 -m src.ga4gh.vrs.extras.vcf_annotation --vcf_in tests/extras/data/test_vcf_input.vcf --vcf_out /data/output.vcf.gz --vrs_file /data/vrs_objects.pkl --seqrepo_root_dir /usr/local/share/seqrepo/2021-01-29
```

To use SeqRepo REST (uses test data for example):

```commandline
python3 -m src.ga4gh.vrs.extras.vcf_annotation --vcf_in tests/extras/data/test_vcf_input.vcf --vcf_out /data/output.vcf.gz --vrs_file /data/vrs_objects.pkl --seqrepo_dp_type rest --seqrepo_base_url http://seqrepo-rest-service:5000/seqrepo
```
