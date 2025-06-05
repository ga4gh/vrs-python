# VRS-in-a-Box VCF Annotator for Assembly GRCh37 or GRCh38
VRS-in-a-Box is a single Docker image that is able to annotate a VCF with VRS IDs
using the `vrs-annotate` tool that is part of `vrs-python`.  The Docker image includes
all the dependencies required for VRS ID computation of genomic variants on one of the
assembled chromosomes for a specific reference assembly.

The Docker image is kept to a minimum by creating an instance of `seqrepo` that
contains only the assembled chromosomes for a single reference assembly.

Prebuilt images are available in Docker Hub: https://hub.docker.com/u/ga4gh

## Using VRS-in-a-Box in Terra
VRS-in-a-Box can be easily added to a workflow in Terra to annotate a VCF file with
VRS IDs.  The `VrsVcfAnnotator.wdl` file contains a simple workflow with one task
that will annotate a VCF file using the pre-built images in Docker Hub.

## Building VRS-in-a-Box
The following instructions describe how to build a VRS-in-a-Box image from scratch.


#### Create the Build Environment
Install any [prerequisites](https://github.com/biocommons/biocommons.seqrepo#requirements)
needed for `seqrepo` and create a Python virtual environment.
```bash
python -m venv venv
source venv/bin/activate
pip install biocommons.seqrepo
```

#### Download the Reference Assembly Sequences
```bash
# GRCh38
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
# GRCh37
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.fna.gz
```

#### Build Seqrepo
```bash
# GRCh38
bash seqrepo.sh GCF_000001405.26_GRCh38_genomic.fna.gz GRCh38
# GRCh37
bash seqrepo.sh GCF_000001405.26_GRCh37_genomic.fna.gz GRCh37
```

## Build Images for Each Assembly
```shell
# GRCh38
docker build --build-arg ASSEMBLY=GRCh38 --build-arg VRS_PYTHON_VERSION=2.1.2 -t ga4gh/vrs-vcf-annotator-grch38:2.1.2 .
# GRCh37
docker build --build-arg ASSEMBLY=GRCh37 --build-arg VRS_PYTHON_VERSION=2.1.2 -t ga4gh/vrs-vcf-annotator-grch37:2.1.2 .
```

## Running the Image in Docker
Run the image to annotate the VCF file `NA12878.vcf` in the current directory:
```shell
docker run -it --rm -v $(pwd):/input ga4gh/vrs-vcf-annotator-grch38:2.1.2 /input/NA12878.vcf --vcf-out /input/NA12878_with_vrs.vcf
```
Run the image to annotate the VCF file `NA12878.vcf` in the current directory and capture the VRS objects in a separate file:
```shell
docker run -it --rm -v $(pwd):/input ga4gh/vrs-vcf-annotator-grch38:2.1.2 /input/NA12878.vcf --vcf-out /input/NA12878_with_vrs.vcf --ndjson-out /input/vrs-objects.json
```
