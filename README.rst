vr-python
!!!!!!!!!!


Installation
@@@@@@@@@@@@

1. Install vr-python
#####################

The following instructions are for Ubuntu 18.04+ and MacOS.
vr-python is unlikely to work on Windows due to dependencies.

::

  git clone --recurse-submodules https://github.com/ga4gh/vr-python.git
  cd vr-python
  python3.7 -m venv venv/3.7
  source venv/3.7/bin/activate
  pip install --upgrade pip setuptools
  pip install -e .
  pip install -e '.[dev,notebooks]'

(Python 3.5 and 3.6 should also work.)


2. Pull seqrepo data
####################

Sequence data are required to normalize sequences and infer `ga4gh`
sequence identifiers.  The notebooks use `SeqRepo
<https://github.com/biocommons/biocommons.seqrepo>`__.  VR
implementers may use SeqRepo or other data source.

Then, download seqrepo with::

  sudo mkdir /usr/local/share/seqrepo
  sudo chown $USER:$USER /usr/local/share/seqrepo
  seqrepo pull

NOTE: This will download approximately 10GB of sequence data.  See
https://github.com/biocommons/biocommons.seqrepo/ for more
information.



Running the Notebooks
@@@@@@@@@@@@@@@@@@@@@

Once installed as described above, type::

  source venv/3.7/bin/activate
  jupyter notebook --notebook-dir notebooks/



Security
@@@@@@@@

A stand-alone security review has been performed on the specification
itself.  This implementation is offered as-is, and without any
security guarantees. It will need an independent security review
before it can be considered ready for use in security-critical
applications. If you integrate this code into your application it is
AT YOUR OWN RISK AND RESPONSIBILITY to arrange for a security audit.
