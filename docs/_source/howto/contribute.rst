.. _contribute:

Contributing: Build, Test and Submit
====================================

If you want to contribute code to POSYDON, the workflow is: get a development
environment, make and test your changes, and open a pull request.
This page focuses on the *practical* side (building the docs locally, running
the tests, running pre-commit). For the overall contribution policy (code of
conduct, PR etiquette) see the :ref:`contributing page <how-to-contribute>`.

Set up a development environment
--------------------------------

Fork and clone the repository, then install with the development extras:

.. code-block:: bash

   git clone https://github.com/your-username/POSYDON.git
   cd POSYDON
   conda create -n posydon-dev python=3.11
   conda activate posydon-dev
   pip install -e .[dev]

Create a branch for your work before you start editing:

.. code-block:: bash

   git checkout -b the-fix-or-feature

Run the test suite
--------------------------

POSYDON uses ``pytest``; you can run the whole suite or a single file:

.. code-block:: bash

   pytest
   pytest posydon/unit_tests/grids/test_psygrid.py

Use the tests as your safety net while you develop, and run them again before
you push.

Run pre-commit
---------------

The repository uses pre-commit hooks (formatting / linting) that must pass.
Install and run them on your changes:

.. code-block:: bash

   pre-commit install
   pre-commit run --all-files

Fix anything the hooks flag before committing.

Build the documentation locally
----------------------------------

The documentation is Sphinx, in ``docs/_source``, built with the
``posydon-debug`` conda environment (which provides the theme and
``nbsphinx``). From ``docs/``:

.. code-block:: bash

   conda activate posydon-debug
   cd docs
   make html

To see warnings as errors (so you do not miss broken references), run a strict
build:

.. code-block:: bash

   sphinx-build -b html -W -n _source _build/check

and check for dead links with:

.. code-block:: bash

   sphinx-build -b linkcheck _source _build/linkcheck

Add or change documentation in the same PR as the code it describes when it
changes user-visible behavior.

Submit a pull request
-------------------------

When your branch is green (tests pass, pre-commit clean, docs build), push it
and open a pull request against the ``main`` branch, describing what you changed
and why. Maintainers will review; please be available to iterate.

Related
--------

* :ref:`how-to-contribute` -- the contribution guidelines and code of conduct.
