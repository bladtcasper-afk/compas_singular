********************************************************************************
Installation
********************************************************************************

:mod:`compas_singular` requires Python 3.9 or later.


Latest
======

The package is not on PyPI yet. Install the latest version from GitHub with pip.

.. code-block:: bash

    pip install "compas_singular @ git+https://github.com/bladtcasper-afk/compas_singular.git@dev"

Optional extras add ``viewer`` for the 3D viewer used by the examples,
``fd`` for constrained smoothing with :mod:`compas_fd`,
and ``rhino`` for the Rhino 8 scene objects.

.. code-block:: bash

    pip install "compas_singular[viewer] @ git+https://github.com/bladtcasper-afk/compas_singular.git@dev"


Development
===========

To work on the code, clone the repository and install it editable from local source,
with the "viewer" extra to run the examples.

.. code-block:: bash

    git clone https://github.com/bladtcasper-afk/compas_singular.git
    cd compas_singular
    pip install -e ".[viewer]"

A ready-made conda environment (Python 3.12, COMPAS, the viewer and an editable install of this package)
is described in ``environment.yml``.

.. code-block:: bash

    conda env create -f environment.yml
    conda activate compas-singular

To build this documentation, install the "dev" extra and run the ``docs`` task.

.. code-block:: bash

    pip install -e ".[dev]"
    invoke docs

The HTML is written to ``dist/docs``.


Verify
======

.. code-block:: bash

    python -c "import compas_singular; print(compas_singular.__file__)"

If this prints ``None`` instead of a path ending in ``src/compas_singular/__init__.py``,
Python found the bare repository folder and treated it as an empty namespace package.
Run one of the install commands above.


Rhino 8
=======

The Rhino commands are in ``rhino_plugin/commands/``, numbered in workflow order
(``CMD00_start`` to ``CMD09_edit_quad_mesh``).
They run on Rhino 8's CPython 3.9 and need :mod:`compas_singular` installed in Rhino's Python.
Close Rhino, then install the package with Rhino's own interpreter.

.. code-block:: bash

    %USERPROFILE%\.rhinocode\py39-rh8\python.exe -m pip install "compas_singular[rhino,fd] @ git+https://github.com/bladtcasper-afk/compas_singular.git@dev"

The modules under :mod:`compas_singular.rhino` need :mod:`compas_rhino` and are only importable inside Rhino.
See :doc:`tutorial/rhino` for the commands and the workflow they implement.
