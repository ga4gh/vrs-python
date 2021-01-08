vrs-python
!!!!!!!!!!


|ci_rel| | |cov| | |pypi_rel|

vrs-python provides Python language support for the [GA4GH Variation
Representation Specification
(VRS)](https://github.com/ga4gh/vr-spec).

This repository contains several related components:

* **ga4gh.vrs package** Python language support for the spec. 

* **ga4gh.vrs.extras package** Python language support for additional
  functionality, including translating from and to other variant
  formats and a REST service to similar functionality.
  `ga4gh.vrs.extras` requires access to supporting data, as described
  below.

* **Jupyter notebooks** Demonstrations of the functionality of
  `ga4gh.vrs` and `ga4gh.vrs.extras` in the form of easy-to-read
  notebooks.



Installing ga4gh.vrs
!!!!!!!!!!!!!!!!!!!!

Installating with pip
@@@@@@@@@@@@@@@@@@@@@

::

   $ pip install ga4gh.vrs[extras]

The `[extras]` argument tells pip to install packages to fullfill the
dependencies of the `ga4gh.vrs.extras` package.


Installing for development
@@@@@@@@@@@@@@@@@@@@@@@@@@

The following instructions are for Ubuntu 18.04+ and MacOS.
vrs-python is unlikely to work on Windows due to dependencies.

::

   $ git clone --recurse-submodules https://github.com/ga4gh/vrs-python.git
   $ cd vrs-python
   $ make devready

(Python 3.5 and 3.6 should also work.)



Installing Dependencies for ga4gh.vrs.extras
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The `ga4gh.vrs.extras` modules are not part of the VR spec per se.
They are bundled with ga4gh.vrs for development and installation
convenience.  These modules depend directly and indrectly on external
data sources of sequences, transcripts, and genome-transcript
alignments.  This section recommends one way to install the biocommons
tools that provide these data.


::

   $ docker volume create --name=uta_vol
   $ docker volume create --name=seqrepo_vol
   $ docker-compose -f misc/stack/docker-compose.yml up

This should start three containers:
* [seqrepo](https://github.com/biocommons/seqrepo): a non-redundant archive of sequences
* [seqrepo-rest-service](https://github.com/biocommons/seqrepo-rest-service): a REST service on seqrepo (localhost:5000)
* [uta](https://github.com/biocommons/uta): a database of transcripts and alignments (localhost:5432)

The seqrepo container will exit as soon as the data are downloaded.

::

   $ docker ps
   CONTAINER ID        IMAGE                                    //  NAMES
   86e872ab0c69        biocommons/seqrepo-rest-service:latest   //  stack_seqrepo-rest-service_1
   a40576b8cf1f        biocommons/uta:uta_20180821              //  stack_uta_1



Running the Notebooks
!!!!!!!!!!!!!!!!!!!!!

Once installed as described above, type::

  $ source venv/3.7/bin/activate
  $ jupyter notebook --notebook-dir notebooks/


The following jupyter extensions are recommended but not required::

  $ pip install jupyter_contrib_nbextensions
  $ jupyter contrib nbextension install --user
  $ jupyter nbextension enable toc2/main
  


Security Note (from the GA4GH Security Team)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

A stand-alone security review has been performed on the specification
itself.  This implementation is offered as-is, and without any
security guarantees. It will need an independent security review
before it can be considered ready for use in security-critical
applications. If you integrate this code into your application it is
AT YOUR OWN RISK AND RESPONSIBILITY to arrange for a security audit.


.. |pypi_rel| image:: https://badge.fury.io/py/ga4gh.vrs.png
  :target: https://pypi.org/project/ga4gh.vrs
  :align: middle

.. |ci_rel| image:: https://travis-ci.org/ga4gh/vrs-python.svg?branch=master
  :target: https://travis-ci.org/ga4gh/vrs-python
  :align: middle 

.. |cov| image:: https://coveralls.io/repos/github/ga4gh/vrs-python/badge.svg?branch=
  :target: https://coveralls.io/github/ga4gh/vrs-python?branch=
