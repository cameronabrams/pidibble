pidibble
========

.. |pypi| image:: https://img.shields.io/pypi/v/pidibble.svg
   :target: https://pypi.org/project/pidibble/
   :alt: PyPI version

.. |pyversions| image:: https://img.shields.io/pypi/pyversions/pidibble.svg
   :target: https://pypi.org/project/pidibble/
   :alt: Supported Python versions

.. |license| image:: https://img.shields.io/pypi/l/pidibble.svg
   :target: https://github.com/cameronabrams/pidibble/blob/main/LICENSE
   :alt: License

.. |tests| image:: https://github.com/cameronabrams/pidibble/actions/workflows/tests.yaml/badge.svg
   :target: https://github.com/cameronabrams/pidibble/actions/workflows/tests.yaml
   :alt: Tests

.. |downloads| image:: https://static.pepy.tech/badge/pidibble
   :target: https://pepy.tech/projects/pidibble
   :alt: PyPI downloads

.. |docs| image:: https://readthedocs.org/projects/pidibble/badge/?version=latest
   :target: https://pidibble.readthedocs.io/en/latest/
   :alt: Documentation status

|pypi| |pyversions| |license| |tests| |downloads| |docs|

**Pidibble** is a Python package for parsing Protein Data Bank (PDB) files in both legacy PDB and modern PDBx/mmCIF formats.  It conforms to the `most recent standard <https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>`_ (v.3.3 Atomic Coordinate Entry Format, ca. 2011).

Unlike parsers like that found in packages like `BioPython <https://biopython.org/wiki/PDBParser>`_, ``pidibble`` provides meaningfully parsed objects from *all* standard PDB record types, not just ``ATOM`` and ``CONECT`` records.

Once installed, the user has access to the :class:`PDBParser` class in the :mod:`pidibble.pdbparse` module.

Pidibble can fetch coordinate files from the `RCSB PDB <https://www.rcsb.org/>`_ and from the `AlphaFold model database <https://alphafold.ebi.ac.uk/>`_.

Pidibble can handle hexadecimal serial numbers in all record types.

New to pidibble?  Start with the :doc:`quickstart`, then work through the
:doc:`User Guide <guide/index>`.  See :ref:`installation` to install the package.

.. note::

   Pidibble is under active development.

.. note::

   Pidibble is used in `pestifer <https://pypi.org/project/pestifer>`_

Contents
--------

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   guide/index
   API <api/API>
   changelog
