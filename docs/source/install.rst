.. _quick_install:

..
   notes:
   * describe Python3 installation at all? Not sure if any prospective VRS user needs their hand held through this, and I worry that we contribute to `this problem <https://xkcd.com/1987/>`_ if we give explicit directions.
   * do we need to explain how to fetch seqrepo data if we also provide instructions for using docker compose? seemingly this belongs in a separate non-quick installation description
   * Is the "send a single SELECT statement to UTA" instruction necessary? Does it make these instructions look deceptively complex?
   * Is the "try each docker image one at a time" suggestion necessary or helpful? Is a person who is unable to think of doing this going to be able to make any further debugging progress from there?

Installation
============

Prerequisites
-------------

While minimal functions of VRS-Python can be fulfilled with just the base Python package, we recommend the following external dependencies to enable the complete functioning of associated tools and utilities.

#. Python >= 3.8

   * Note: Python 3.10 is required for developers contributing to VRS-Python.

#. `PostgreSQL <PostgreSQL>`_, and ``libpq`` developer files

.. tabs::

   .. tab:: MacOS

      Use Homebrew to install PostgreSQL and ``libpq``.

      .. code-block:: shell

         brew update  # make sure Homebrew is up to date
         brew install libpq
         brew install postgresql@14

      A Python interpreter can be installed from the `official Python website <https://www.python.org/downloads/>`_. Alternatively, it can be installed via Homebrew (see `this writeup <https://realpython.com/installing-python/#how-to-install-python-on-macos>`_ for more information):

      .. code-block:: shell

         brew install python3

      See the `Homebrew documentation <https://docs.brew.sh/Installation>`_ for more information and help with troubleshooting.

   .. tab:: Ubuntu

      Use ``apt`` to install the prerequisites. See the `documentation <https://ubuntu.com/server/docs/package-management>`_ for more information.

      .. code-block:: shell

         sudo apt install gcc libpq-dev python3-dev

Setting up VRS-Python
---------------------

VRS-Python is available on `PyPI <https://pypi.org/project/ga4gh.vrs/>`_.

.. code-block:: shell

   pip install 'ga4gh.vrs[extras]'

Calling  ``pip`` with the ``[extras]`` option installs dependencies needed for use of the ``ga4gh.vrs.extras`` package, such as VCF annotation and variant translation.

Setting up external data sources
--------------------------------

In addition to the reference implementation of VRS provided in ``ga4gh.vrs``, VRS-Python includes additional tools to support core sequence lookup and annotation functions, available in the ``ga4gh.vrs.extras`` package. These modules depend directly and indirectly on external data sources of sequences, transcripts, and genome-transcript alignments.

Fetching SeqRepo data
+++++++++++++++++++++

`SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_ provides fast access to biological sequences and sequence metadata. Install SeqRepo and use the provided command-line tools to fetch a recent data snapshot:

.. code-block:: shell

   pip install seqrepo
   export SEQREPO_VERSION=2021-01-29  # or newer if available -- check `seqrepo list-remote-instances`
   sudo mkdir /usr/local/share/seqrepo
   sudo chown $USER /usr/local/share/seqrepo
   seqrepo pull -i $SEQREPO_VERSION

.. note::
   Users may encounter permission errors like the following while executing ``seqrepo pull``: ::

      PermissionError: [Error 13] Permission denied: '/usr/local/share/seqrepo/2021-01-29._fkuefgd' -> '/usr/local/share/seqrepo/2021-01-29'

   If so, use ``sudo`` to move the data manually: ::

      sudo mv /usr/local/share/seqrepo/$SEQREPO_VERSION.* /usr/local/share/seqrepo/$SEQREPO_VERSION

Launching remaining services with Docker
++++++++++++++++++++++++++++++++++++++++

Installation of the remaining dependencies, `SeqRepo REST Service <https://github.com/biocommons/seqrepo-rest-service>`_ and the `Universal Transcript Archive <https://github.com/biocommons/uta>`_ is simplest by way of their provided `Docker <https://www.docker.com/>`_ images.

.. note::
   If you would like to use local instances of UTA, see `UTA <https://github.com/biocommons/uta>`_ directly. We do provide some additional setup help [here](./docs/setup_help/).

Launch the requisite Docker services. Depending on your network and hoost, the *first* run is likely to take several minutes in order to fully download and install all data. Subsequent startups should be nearly instantaneous.

.. code-block:: shell

   docker volume create --name=uta_vol
   docker volume create --name=seqrepo_vol
   docker-compose up

After all images are acquired, verify that the SeqRepo REST Service and UTA containers are running by calling the ``docker ps`` command.

.. code-block:: shell

   $ docker ps
   CONTAINER ID        IMAGE                                    //  NAMES
   86e872ab0c69        biocommons/seqrepo-rest-service:latest   //  vrs-python_seqrepo-rest-service_1
   a40576b8cf1f        biocommons/uta:uta_20210129b              //  vrs-python_uta_1

You can test UTA and seqrepo installations like so: ::

   $ psql -XAt postgres://anonymous@localhost/uta -c 'select count(*) from transcript'
   # the result should be somewhere in the vicinity of 24909


It doesn't work
_______________

Here are some things to try.

- Bring up one service at a time. For example, if you haven't download seqrepo
  yet, you might see this: ::

   $ docker-compose up seqrepo-rest-service
   Starting vrs-python_seqrepo-rest-service_1 ... done
   Attaching to vrs-python_seqrepo-rest-service_1
   seqrepo-rest-service_1  | 2022-07-26 15:59:59 seqrepo_rest_service.__main__[1] INFO Using seqrepo_dir='/usr/local/share/seqrepo/2021-01-29' from command line
   ⋮
   seqrepo-rest-service_1  | OSError: Unable to open SeqRepo directory /usr/local/share/seqrepo/2021-01-29
   vrs-python_seqrepo-rest-service_1 exited with code 1
