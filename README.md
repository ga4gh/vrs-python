# vrs-python

**VRS-Python** provides Python language support and a reference implementation for the
[GA4GH Variation Representation Specification(VRS)](https://github.com/ga4gh/vrs).

## Information

[![license](https://img.shields.io/badge/license-Apache-green)](https://github.com/ga4gh/vrs-python/blob/main/LICENSE) [![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ga4gh/vrs-python/main?labpath=notebooks)

## Releases

[![gitHub tag](https://img.shields.io/github/v/tag/ga4gh/vrs-python.svg)](https://github.com/ga4gh/vrs-python/releases) [![pypi](https://img.shields.io/pypi/v/ga4gh.vrs.svg)](https://pypi.org/project/ga4gh.vrs/)

## Development

 [![action status](https://github.com/ga4gh/vrs-python/actions/workflows/python-cqa.yml/badge.svg)](https://github.com/ga4gh/vrs-python/actions/workflows/python-cqa.yml) [![issues](https://img.shields.io/github/issues-raw/ga4gh/vrs-python.svg)](https://github.com/ga4gh/vrs-python/issues)
[![GitHub Open Pull Requests](https://img.shields.io/github/issues-pr/ga4gh/vrs-python.svg)](https://github.com/ga4gh/vrs-python/pull/) [![GitHub license](https://img.shields.io/github/contributors/ga4gh/vrs-python.svg)](https://github.com/ga4gh/vrs-python/graphs/contributors/) [![GitHub stars](https://img.shields.io/github/stars/ga4gh/vrs-python.svg?style=social&label=Stars)](https://github.com/ga4gh/vrs-python/stargazers) [![GitHub forks](https://img.shields.io/github/forks/ga4gh/vrs-python.svg?style=social&label=Forks)](https://github.com/ga4gh/vrs-python/network)

## Features

- Pydantic implementation of GKS core models and VRS models
- Algorithm for generating consistent, globally unique identifiers for variation without a central authority
- Algorithm for performing fully justified allele normalization
- Translating from and to other variant formats
- [Annotate VCFs with VRS](https://github.com/ga4gh/vrs-python/blob/main/docs/extras/vcf_annotator.md)
- Convert GA4GH objects between inlined and referenced forms

## Known Issues

**You are encouraged to** [browse issues](https://github.com/ga4gh/vrs-python/issues).
All known issues are listed there. Please report any issues you find.

## Installing VRS-Python Locally

### Prerequisites

- Python >= 3.10
  - _Note: Python 3.12 is required for developers contributing to VRS-Python_
- libpq
- postgresql

#### MacOS

You can use Homebrew to install the prerequisites. See the
[Homebrew documentation](https://docs.brew.sh/Installation) for how to install. Make
 sure Homebrew is up-to-date by running `brew update`.

```shell
brew install libpq
brew install python3
brew install postgresql@14
```

#### Ubuntu

```shell
sudo apt install gcc libpq-dev python3-dev
```

### Installation Steps

#### 1. Install VRS-Python with pip

VRS-Python is available on [PyPI](https://pypi.org/project/ga4gh.vrs/).

```shell
pip install 'ga4gh.vrs[extras]'
```

The `[extras]` argument tells pip to install packages to fulfill the dependencies of the
`ga4gh.vrs.extras` package.

#### 2. Install External Data Sources

The `ga4gh.vrs.extras` modules are not part of the VR spec per se. They are
bundled with ga4gh.vrs for development and installation convenience. These
modules depend directly and indirectly on external data sources of sequences,
transcripts, and genome-transcript alignments.

First, you must install a local [SeqRepo](https://github.com/biocommons/biocommons.seqrepo):

```shell
pip install seqrepo
export SEQREPO_VERSION=2024-02-20  # or newer if available -- check `seqrepo list-remote-instances`
sudo mkdir -p /usr/local/share/seqrepo
sudo chown $USER /usr/local/share/seqrepo
seqrepo pull -i $SEQREPO_VERSION
seqrepo update-latest
```

If you encounter a permission error similar to the one below:

```shell
PermissionError: [Error 13] Permission denied: '/usr/local/share/seqrepo/2024-02-20._fkuefgd' -> '/usr/local/share/seqrepo/2024-02-20'
```

Try moving data manually with `sudo`:

```shell
sudo mv /usr/local/share/seqrepo/$SEQREPO_VERSION.* /usr/local/share/seqrepo/$SEQREPO_VERSION
```

To make installation easy, we recommend using Docker to install the other Biocommons
tools - [SeqRepo REST](https://github.com/biocommons/seqrepo-rest-service) and
[UTA](https://github.com/biocommons/uta). If you would like to use local instances of UTA,
see [UTA](https://github.com/biocommons/uta) directly. We do provide some additional
setup help [here](./docs/setup_help/).

Next, run the following commands:

```shell
docker volume create --name=uta_vol
docker volume create --name=seqrepo_vol
docker-compose up
```

This should start three containers:

- [seqrepo](https://github.com/biocommons/seqrepo): downloads seqrepo into a
  docker volume and exits
- [seqrepo-rest-service](https://github.com/biocommons/seqrepo-rest-service): a
  REST service on seqrepo (localhost:5000)
- [uta](https://github.com/biocommons/uta): a database of transcripts and
  alignments (localhost:5432)

Check that the containers are running, by running:

```shell
$ docker ps
CONTAINER ID        IMAGE                                    //  NAMES
86e872ab0c69        biocommons/seqrepo-rest-service:latest   //  vrs-python_seqrepo-rest-service_1
a40576b8cf1f        biocommons/uta:uta_20210129b              //  vrs-python_uta_1
```

Depending on your network and host, the _first_ run is likely to take 5-15
minutes in order to download and install data. Subsequent startups should be
nearly instantaneous.

You can test UTA and seqrepo installations like so:

```shell
$ psql -XAt postgres://anonymous@localhost/uta -c 'select count(*) from uta_20210129b.transcript'
314227
```

##### It doesn't work

Here are some things to try.

- Bring up one service at a time. For example, if you haven't download seqrepo
  yet, you might see this:

  ```shell
  $ docker-compose up seqrepo-rest-service
  Starting vrs-python_seqrepo-rest-service_1 ... done
  Attaching to vrs-python_seqrepo-rest-service_1
  seqrepo-rest-service_1  | 2022-07-26 15:59:59 seqrepo_rest_service.__main__[1] INFO Using seqrepo_dir='/usr/local/share/seqrepo/2024-02-20' from command line
  ⋮
  seqrepo-rest-service_1  | OSError: Unable to open SeqRepo directory /usr/local/share/seqrepo/2024-02-20
  vrs-python_seqrepo-rest-service_1 exited with code 1
  ```

## VRS-Python and VRS Version Correspondence

The ga4gh/vrs-python repo embeds the ga4gh/vrs repo as a git submodule for testing purposes.
Each ga4gh.vrs package on PyPI embeds a particular version of VRS. The
correspondences between the packages that are **currently maintained** may be summarized as:

| vrs-python branch | vrs-python tag/version | vrs branch | vrs version |
| --- | --- | --- | --- |
| [main](https://github.com/ga4gh/vrs-python/tree/main) _(default branch)_ | 2.x | [2.x](https://github.com/ga4gh/vrs/tree/2.x) | 2.x |
| [1.x](https://github.com/ga4gh/vrs-python/tree/1.x) | 0.8.x | [1.x](https://github.com/ga4gh/vrs/tree/1.x) | 1.x |

⚠ **Note: Only 2.x branch is being actively maintained. The 1.x branch will only be maintained for bug fixes.**

⚠ **Developers: See the development section below for recommendations for using submodules
gracefully (and without causing problems for others!).**

### Previous VRS-Python and VRS Version Correspondence

The correspondences between the packages that are **no longer maintained** may be summarized as:

| vrs-python branch | vrs-python tag/version | vrs branch | vrs version |
| --- | --- | --- | --- |
| [0.9](https://github.com/ga4gh/vrs-python/tree/0.9) | 0.9.x | [metaschema-update](https://github.com/ga4gh/vrs/tree/metaschema-update) | N/A |
| [0.7](https://github.com/ga4gh/vrs-python/tree/0.7) | 0.7.x | [1.2](https://github.com/ga4gh/vrs/tree/1.2) | 1.2.x |
| [0.6](https://github.com/ga4gh/vrs-python/tree/0.6) | 0.6.x | [1.1](https://github.com/ga4gh/vrs/tree/1.1) | 1.1.x |

## Developers

This section is intended for developers who contribute to VRS-Python.

### Installing for development

Fork the repo at <https://github.com/ga4gh/vrs-python/>.

```shell
git clone --recurse-submodules git@github.com:YOUR_GITHUB_ID/vrs-python.git
cd vrs-python
make devready
source venv/3.12/bin/activate
```

If you already cloned the repo, but forgot to include `--recurse-submodules` you can run:

```shell
git submodule update --init --recursive
```

#### Submodules

vrs-python embeds vrs as a submodule, only for testing purposes. When checking out vrs-python and switching
branches, it is important to make sure that the submodule tracks vrs-python
correctly. The recommended way to do this is `git config --global submodule.recurse true`. **If you don't set submodule.recurse, developers and
reviewers must be extremely careful to not accidentally upgrade or downgrade
schemas with respect to vrs-python.**

Alternatively, see `misc/githooks/`.

### Testing

This package implements typical unit tests for ga4gh.core and ga4gh.vrs. This
package also implements the compliance tests from vrs (vrs/validation) in the
tests/validation/ directory.

To run tests:

```shell
make test
```

## Running the Notebooks

The notebooks **do not** require you to setup SeqRepo or UTA from
[Install External Data Sources](#2-install-external-data-sources).

### Running the Notebooks on Binder

[Binder](https://mybinder.org/) allows you to create custom computing environments that can be shared and used by many remote users.

You can access the notebooks on Binder [here](https://mybinder.org/v2/gh/ga4gh/vrs-python/main?labpath=notebooks).

### Running the Notebooks on the Terra platform

[Terra](https://terra.bio) is a cloud platform for biomedical research developed by the Broad Institute, Microsoft and Verily. The platform includes preconfigured environments that provide user-friendly access to various applications commonly used in bioinformatics, including Jupyter Notebooks.

We have created a public [`VRS-demo-notebooks`](https://app.terra.bio/#workspaces/terra-outreach/VRS-demo-notebooks) workspace in Terra that contains the demo notebooks along with instructions for running them with minimal setup. To get started, see either the [`VRS-demo-notebooks`](https://app.terra.bio/#workspaces/terra-outreach/VRS-demo-notebooks) workspace or the [`Terra.ipynb`](notebooks/Terra.ipynb) notebook in this repository.

### Running the Notebooks with VS Code

[VS Code](https://code.visualstudio.com/) is a code editor developed by Microsoft. It is lightweight, highly customizable, and supports a wide range of programming languages, with a robust extension system. You can download VS Code [here](https://code.visualstudio.com/Download).

1. Open VS Code.
2. Use Extensions view (Ctrl+Shift+X or ⌘+Shift+X) to install the [Jupyter extension](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter).
3. Navigate to your vrs-python project folder and open it in VS Code.
4. In a notebook, click `Select Kernel` at the top right. Select the option where the path is `venv/3.12/bin/python3`. See [here](https://code.visualstudio.com/docs/datascience/jupyter-kernel-management) for more information on managing Jupyter Kernels in VS Code.
5. After selecting the kernel you can now run the notebook.

## Security Note (from the GA4GH Security Team)

A stand-alone security review has been performed on the specification itself.
This implementation is offered as-is, and without any security guarantees. It
will need an independent security review before it can be considered ready for
use in security-critical applications. If you integrate this code into your
application it is AT YOUR OWN RISK AND RESPONSIBILITY to arrange for a security
audit.
