.. _installation:

Installation
============

To use ``Pidibble``, install it from PyPI:

.. code-block:: bash

   (.venv) $ pip install pidibble

Reading PDBx/mmCIF files needs one more package, the wwPDB ``mmcif`` reader.
It is an optional extra because it is a compiled package that not every
channel carries; PDB parsing, writing and citations work without it:

.. code-block:: bash

   (.venv) $ pip install 'pidibble[mmcif]'

Without it, ``PDBParser(..., input_format='mmCIF').parse()`` raises an
``ImportError`` naming the fix, before anything is downloaded.

Pidibble is under active development, so you can also install the latest version from the GitHub repository:

.. code-block:: bash

   (.venv) $ pip install git+https://github.com/cameronabrams/pidibble.git