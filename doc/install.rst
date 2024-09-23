
===============================================================================
Installation
===============================================================================

There are different ways to install `cocoatree`:

- :ref:`Install the latest official release. <install_official_release>` This
  is the best approach for most users and will provide a stable version.
- :ref:`Building the package from source. <install_from_source>` This is best
  for users who want the latest-and-greatest features and aren't afraid of
  running brand-new code. This is also needed for users who wish to contribute
  to the project.

.. _install_official_release:

Install the latest release
==========================

Install a version of Python 3, for instance from https://www.python.org.

Then run::

  pip install -U cocoatree

In order to check your installation you can use::

  python -m pip show cocoatree # to see which version and where cocoatree is installed
  python -m pip freeze # to see all packages installed in the active environment.


Note that in order to avoid potential conflicts with other packages it is
strongly recommended to use a virtual environment, e.g. python3 `virtualenv`
or conda environments.

.. _install_from_source:


Building the package from source
================================

Building from source is required to work on a contribution (bug fix, new
feature, code or documentation improvement).

1. Use Git to check out the latest source from the cocoatree repository on Github.::

    git clone git://github.com/tree-timc/cocoatree.git
    cd cocoatree

   If you plan on submitting a pull-request, you should clone from your fork instead.

2. Install the following dependencies if they are not already installed:

   - Python (>= 3.11)
   - NumPy (>= 2.1),
   - SciPy (>= 0.1dd9.1),
   - Matplotlib (>= 3.9),
   - biopython

2. Optional (but recommended): create and activate a dedicated virtualenv or conda environment.

3. Build the project with pip in Editable mode::

    pip install --verbose --editable .

   If you run the development version, it is cumbersome to reinstall the
   package each time you update the sources. Therefore it is recommended that
   you install in with the `pip install --editable .` command, which allows
   you to edit the code in-place. This builds the extension in place and
   creates a link to the development directory (see the pip docs).

   This is fundamentally similar to using the command `python setup.py develop`
   (see the setuptool docs). It is however preferred to use pip.

   On Unix-like systems, you can equivalently type make in from the top-level
   folder. Have a look at the Makefile for additional utilities.

4. Check that the installed `cocoatree` has a version number ending with `.dev0`::

      python -c "import cocoatree; cocoatree.__version__"


Additional dependencies for building the documentation
------------------------------------------------------

In addition to the formentioned dependencies for `cocoatree`, in order to
build the documentation you will need:

- sphinx (>= 1.0)
- sphinx-gallery
- numpydoc
- pillow

These can be installed with pip using the command::

  pip install -U sphinx sphinx-gallery numpydoc

Additional dependencies for running the tests
---------------------------------------------

In addition to the formentioned dependencies for `cocoatree`, in order to
build the documentation you will need:

- pytest
- pytest-cov
- flake8

These can be installed with pip using the following command::

  pip install -U pytest pytest-cov


