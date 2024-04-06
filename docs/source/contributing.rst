.. _contributing:

Contributing
============

Bug reports and feature requests
--------------------------------

Bugs and new feature requests can be submitted to the `issue tracker on GitHub <https://github.com/ga4gh/vrs-python/issues>`_. See `this StackOverflow post <https://stackoverflow.com/help/minimal-reproducible-example>`_ for tips on how to craft a helpful bug report.

Development setup
-----------------

Clone the repository: ::

    git clone https://github.com/ga4gh/vrs-python
    cd vrs-python

Then initialize a virtual environment: ::

    python3 -m virtualenv venv
    source venv/bin/activate
    python3 -m pip install -e '.[dev,extras,notebooks,docs]'

Tests
-----

Tests are executed with `pytest <https://docs.pytest.org/en/7.1.x/getting-started.html>`_ via a Makefile target: ::

    make test

Documentation
-------------

The documentation is built with Sphinx, which is included as part of the ``docs`` dependency group. Navigate to the ``docs/`` subdirectory and use ``make`` to build the HTML version: ::

    cd docs
    make html

See the `Sphinx documentation <https://www.sphinx-doc.org/en/master/>`_ for more information.
