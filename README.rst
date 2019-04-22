vmc-python
!!!!!!!!!!

VMC is the Variation Modelling Collaboration.  The VMC's mission is to
standardize the exchange of biological sequence variants among
computer systems.  The primary "products" of this effort are:

  #. A terminology document defining core elements of the data model.

  #. Machine-readable specifications for the data model.

This repository contains code to demonstrate the use fo the VMC data
model.


**NOTE:** This project is in-progress.  


Installation
@@@@@@@@@@@@

1. Install vmc-python
#####################

The following instructions are for Ubuntu 18.04+ and MacOS.
vmc-python is unlikely to work on Windows due to dependencies.

::

  git clone --recurse-submodules https://github.com/ga4gh/vmc-python.git
  cd vmc-python
  python3.7 -m venv venv/3.7
  source venv/3.7/bin/activate
  pip install --upgrade pip setuptools
  pip install -e .
  pip install -e '.[dev,notebooks]'

(Python 3.5 and 3.6 should also work.)


2. Pull seqrepo data
####################

Sequence data are required to normalize sequences and infer VMC
sequence identifiers.  The notebooks use `SeqRepo
<https://github.com/biocommons/biocommons.seqrepo>`__.  VMC
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
