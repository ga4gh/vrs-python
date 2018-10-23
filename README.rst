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


Installation (required to use notebooks)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Ubuntu 18.04+
###############

::

  git clone --recurse-submodules https://github.com/ga4gh/vmc-python.git
  make venv/3.6
  source venv/3.6/bin/activate
  make install


MacOS
########

(Please contribute instructions)


Windows
#######

You need a different kind of help.


Running the Notebooks
@@@@@@@@@@@@@@@@@@@@@

Once installed as described above, type::

  source venv/3.6/bin/activate
  pip install jupyter
  jupyter notebook --notebook-dir notebooks/
