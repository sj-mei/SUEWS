.. _installation:


Installation
============



Formal releases
---------------

Since 2023, SUEWS is available as a command line tool via its Python wrapper package `SuPy (SUEWS in Python) <SuPy>`_ on `PyPI`_ and `conda-forge`_.

Installing Python
*****************
These instructions will set you up with `miniforge <https://conda-forge.org/docs/user/introduction.html>`__, which makes it easy to install and manage Python packages.

To install the miniforge Python distribution follow `the miniforge installation instructions <https://github.com/conda-forge/miniforge#install>`__.

This makes installing ``supy`` and many other packages in the scientific Python ecosystem much easier and quicker.
It also provides many pre-compiled binaries that are not available on PyPI.

Installing SuPy
***************

One can install SUEWS using ``pip``:

.. code-block:: bash

    pip install supy

or ``conda``:

.. code-block:: bash

    conda install -c conda-forge supy


.. tip::

    ``conda`` could be slower than ``pip`` to install SUEWS.
    If you are using ``conda``, please consider using ``mamba`` instead:

    .. code-block:: shell

        mamba install -c conda-forge supy

    ``mamba`` is a drop-in replacement for ``conda`` that is faster and solves some common problems with ``conda``.
    More details about ``mamba`` can be found at `mamba`_.



.. note::

    The Fortran-based binaries build prior to 2023 are still available at the `SUEWS download page`_.
    However, they are not maintained anymore so users are encouraged to use the Python-based packages instead.



.. _PyPI: https://pypi.org/project/supy/
.. _conda-forge: https://anaconda.org/conda-forge/supy
.. _mamba: https://github.com/mamba-org/mamba
.. _SuPy: :ref:`supy_index`



Development build
-----------------

.. warning::

The development build can be highly unstable and is not recommended for production use.
However, it is automatically constructed every week for testing purposes and we are happy to receive feedback on the development build.


To install the development build of SUEWS, you need to install ``supy`` in the development mode:

1. git clone the repository::

    git clone https://github.com/UMEP-dev/SUEWS.git

2. navigate to the directory of the cloned repository::

    cd SUEWS

3. install the package in the development mode::

    pip install -e src/supy


