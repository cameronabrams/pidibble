Welcome to Pidibble's documentation!
====================================

**Pidibble** is a Python package for parsing Protein Data Bank (PDB) files.  It conforms to the `most recent standard <https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>`_ (v.3.3 Atomic Coordinate Entry Format, ca. 2011).

Unlike parsers like that found in packages like `BioPython <https://biopython.org/wiki/PDBParser>`_, ``pidibble`` provides meaningfully parsed objects from *all* standard PDB record types, not just ATOMs and CONECTs.

Once installed, the user has access to the ``PDBParser`` class in the ``pidibble.pdbparser`` module.

Pidibble can fetch coordinate files from the RSCB PDB and from the AlphaFold model database.

Pidibble can handle hexadecimal serial numbers in all record types.

Check out the :doc:`usage` section for further information, including on :ref:`installation` of the package.

.. note::

   Pidibble is under active development.

.. note::

   Pidibble is used in `pestifer <https://pypi.org/project/pestifer>`_

Contents
--------

.. toctree::

   usage
   api
