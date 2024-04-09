.. _install:

Installation
============

VRS-Python can be installed from `PyPI <https://pypi.org/projects/ga4gh.vrs>`_. :ref:`Quick install <quick-install>` instructions are available for users who are able to launch `Docker <https://www.docker.com/>`_ instances in their local environment. Otherwise, the :ref:`full install <full-install>` instructions describe how to stand up local copies of external data dependencies (`SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_ and the `Universal Transcript Archive <https://github.com/biocommons/uta>`_) and make them available to VRS-Python.

.. note::

   The core models and functions of VRS-Python do not require any external data sources. Users who don't plan on employing any of the annotation and translation functions provided in the `extras` module can simply install VRS-Python with the base dependency set:

   .. code-block:: shell

      $ python3 -m pip install ga4gh.vrs

Prerequisites
-------------

#. Python >= 3.8

   * Note: Python 3.10 is required for developers contributing to VRS-Python.

#. `PostgreSQL <PostgreSQL>`_, and ``libpq`` developer files:

.. tabs::

   .. tab:: MacOS

      Use Homebrew to install PostgreSQL and ``libpq``.

      .. code-block:: shell

         brew update  # make sure Homebrew is up to date
         brew install libpq
         brew install postgresql@14

      See the `Homebrew documentation <https://docs.brew.sh/Installation>`_ for more information and help with troubleshooting.

   .. tab:: Ubuntu

      Use ``apt`` to install the prerequisites. See the `documentation <https://ubuntu.com/server/docs/package-management>`_ for more information.

      .. code-block:: shell

         sudo apt install gcc libpq-dev python3-dev

.. _quick-install:

Quick install
-------------

First, install VRS-Python and the ``extras`` dependency group via `PyPI`_.

.. code-block:: shell

   pip install 'ga4gh.vrs[extras]'

UTA and SeqRepo can be launched for use in VRS-Python by way of the `Docker Compose <https://docs.docker.com/compose/>`_ file provided in the VRS-Python source repository. First, acquire the file:

.. code-block:: shell

   curl -O https://github.com/ga4gh/vrs-python/blob/main/docker-compose.yml

Then, launch the requisite Docker services. Depending on your network and host, the *first* run is likely to take several minutes in order to fully download and install all data. Subsequent startups should be nearly instantaneous.

.. code-block:: shell

   docker volume create --name=uta_vol
   docker volume create --name=seqrepo_vol
   docker-compose up

After all images are acquired, verify that the SeqRepo REST Service and UTA containers are running by calling the ``docker ps`` command.

.. code-block:: shell

   $ docker ps
   CONTAINER ID        IMAGE                                    //  NAMES
   86e872ab0c69        biocommons/seqrepo-rest-service:latest   //  vrs-python_seqrepo-rest-service_1
   a40576b8cf1f        biocommons/uta:uta_20210129b             //  vrs-python_uta_1

.. code-block:: shell

   pip install 'ga4gh.vrs[extras]'


.. _full-install:

Full install
------------

Users unable to run Docker on their local machines, or intending to acquire direct access to the sequence and transcript data utilized in VRS-Python's translation tools, may follow these steps

First, install VRS-Python and the ``extras`` dependency group via `PyPI`_.

.. code-block:: shell

   pip install 'ga4gh.vrs[extras]'

Next, set up and acquire a recent SeqRepo data snapshot.

.. code-block:: shell

   export SEQREPO_VERSION=2021-01-29  # or newer, if available -- check with the `seqrepo list-remote-instances` command
   sudo mkdir /usr/local/share/seqrepo
   sudo chown $USER /usr/local/share/seqrepo
   seqrepo pull -i $SEQREPO_VERSION

Users frequently report encountering a `PermissionError` while calling the `pull` command. If necessary, manually move data from its temporary file location:

.. code-block:: shell

   sudo mv /usr/local/share/seqrepo/$SEQREPO_VERSION.* /usr/local/share/seqrepo/$SEQREPO_VERSION

Next, install UTA from a database dump. As noted in the `UTA documentation <https://github.com/biocommons/uta?tab=readme-ov-file#installing-from-database-dumps>`_, PostgreSQL data and binary locations can vary substantially between system architectures, operating systems, so users may need to adapt these instructions to their own working environments as needed.

To begin, start PostgreSQL service:

.. tabs::

   .. tab:: ARM64 MacOS

      .. code-block:: shell

         pg_ctl -D /opt/homebrew/var/postgres start

   .. tab:: Intel MacOS

      .. code-block:: shell

         pg_ctl -D /usr/local/var/postgresql@14 start

   .. tab:: Ubuntu

      TODO

Next, create users and the database needed for UTA database access:

.. code-block::

   createuser -U postgres uta_admin
   createuser -U postgres anonymous
   createdb -U postgres -O uta_admin uta

Then acquire and load the most recent UTA database dump:

.. code-block::

   uta_v=uta_20210129b
   gzip -cdq $uta_v.pgd.gz | psql -U uta_admin -1 -v ON_ERROR_STOP=1 -d uta -Eae
