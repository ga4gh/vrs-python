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

Ubuntu 18.04+
$$$$$$$$$$$$$

::

  git clone --recurse-submodules https://github.com/ga4gh/vmc-python.git
  python3.6 -m venv venv/3.6
  source venv/3.6/bin/activate
  pip install --upgrade pip setuptools
  pip install -e .
  pip install -e '.[notebooks]'


MacOS
$$$$$

::

  git clone --recurse-submodules https://github.com/ga4gh/vmc-python.git
  python3.7 -m venv vmc-python
  cd vmc-python
  pip3 install --upgrade pip setuptools
  pip3 install -e .
  pip3 install -e '.[notebooks]'


Windows
$$$$$$$

You need a different kind of help.



  seqrepo pull


2. Pull seqrepo data
####################

Sequence data are required to normalize sequences and infer VMC
sequence identifiers.  The notebooks use `SeqRepo
<https://github.com/biocommons/biocommons.seqrepo>`__.  VMC
implementers may use SeqRepo or other data source.


Running the Notebooks
@@@@@@@@@@@@@@@@@@@@@

Once installed as described above, type::

  source venv/3.6/bin/activate
  jupyter notebook --notebook-dir notebooks/
