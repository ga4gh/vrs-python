# vr-python

vr-python provides Python language support for the [GA4GH Variation
Representation Specification
(vr-spec)](https://github.com/ga4gh/vr-spec).

This repository contains several related components:

* **ga4gh.vr package** Python language support for the spec. 

* **ga4gh.vr.extras package** Python language support for additional
  functionality, including translating from and to other variant
  formats and a REST service to similar functionality.
  `ga4gh.vr.extras` requires access to supporting data, as described
  below.

* **Jupyter notebooks** Demonstrations of the functionality of
  `ga4gh.vr` and `ga4gh.vr.extras` in the form of easy-to-read
  notebooks.



# Installing ga4gh.vr

## Using pip

```
pip install ga4gh.vr
```

## Installing for development

The following instructions are for Ubuntu 18.04+ and MacOS.
vr-python is unlikely to work on Windows due to dependencies.

```
git clone --recurse-submodules https://github.com/ga4gh/vr-python.git
cd vr-python
make devready
```

(Python 3.5 and 3.6 should also work.)



# Installing Dependencies for ga4gh.vr.extras

The `ga4gh.vr.extras` modules are not part of the VR spec per se.
They are bundled with ga4gh.vr for development and installation
convenience.  These modules depend directly and indrectly on external
data sources of sequences, transcripts, and genome-transcript
alignments.  This section recommends one way to install the biocommons
tools that provide these data.


```
docker-compose -f misc/stack/docker-compose.yml up
```

This should start three containers:
* [seqrepo](https://github.com/biocommons/seqrepo): a non-redundant archive of sequences
* [seqrepo-rest-service](https://github.com/biocommons/seqrepo-rest-service): a REST service on seqrepo (localhost:5000)
* [uta](https://github.com/biocommons/uta): a database of transcripts and alignments (localhost:5432)



# Running the Notebooks

Once installed as described above, type:

```
source venv/3.7/bin/activate
jupyter notebook --notebook-dir notebooks/
```


# Security Note (from the GA4GH Security Team)

A stand-alone security review has been performed on the specification
itself.  This implementation is offered as-is, and without any
security guarantees. It will need an independent security review
before it can be considered ready for use in security-critical
applications. If you integrate this code into your application it is
AT YOUR OWN RISK AND RESPONSIBILITY to arrange for a security audit.
