Welcome to Pidibble's documentation!
===================================

**Pidibble** is a Python package for parsing standard Protein Data Bank (PDB) files.  It conforms to the `most recent standard <https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>`_ (v.3.3 Atomic Coordinate Entry Format, ca. 2011).

Unlike parsers like that found in packages like `BioPython <https://biopython.org/wiki/PDBParser>`_, ``pidibble`` provides meaningfully parsed objects from *all* standard PDB record types, not just ATOMs and CONECTs.

Once installed, the user has access to the ``PDBParser`` class in the ``pidibble.pdbparser`` module.

Check out the :doc:`usage` section for further information, including
on :ref:`installation` of the packge.

.. note::

   Pidibble is under active development.

Contents
--------

.. toctree::

   usage
   api
