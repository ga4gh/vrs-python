.. _install:

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

      You can use Homebrew to install the prerequisites. See the `Homebrew documentation <https://docs.brew.sh/Installation>`_ for how to install.

      .. code-block:: shell

         brew update  # make sure Homebrew is up to date
         brew install libpq
         brew install python3
         brew install postgresql@14

   .. tab:: Ubuntu

      .. code-block:: shell

         sudo apt install gcc libpq-dev python3-dev

Setting up VRS-Python
---------------------

VRS-Python is available on `PyPI <https://pypi.org/project/ga4gh.vrs/>`_.

.. code-block:: shell

   pip install 'ga4gh.vrs[extras]'

The ``[extras]`` option tells ``pip`` to install packages to fulfill the dependencies of the
``ga4gh.vrs.extras`` package.

Setting up external data sources
--------------------------------

In addition to the reference implementation of VRS provided in ``ga4gh.vrs``, VRS-Python includes additional tools to support core sequence lookup and annotation functions, available in the ``ga4gh.vrs.extras`` package. These modules depend directly and indirectly on external data sources of sequences, transcripts, and genome-transcript alignments.

Fetching SeqRepo data
+++++++++++++++++++++

.. warning::
   TODO is this section necessary? could we just recommend docker?

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
   TODO fill in info about local data usage.... If you would like to use local instances of UTA, see [UTA](https://github.com/biocommons/uta) directly. We do provide some additional setup help [here](./docs/setup_help/).

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

.. warning::
   the rest of this stuff doesn't actually seem super helpful or necessary.

   - the "UTA and seqrepo installation test" only tests UTA

   - The "try each one at a time" tip is fine I guess but is probably not going to help anyone who wasn't able to think of that themselves

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
