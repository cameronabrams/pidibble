Loading structures
==================

Every use of pidibble starts by constructing a
:class:`~pidibble.pdbparse.PDBParser` and calling
:meth:`~pidibble.pdbparse.PDBParser.parse`.  The constructor decides *where* the
data comes from and *what format* it is in; ``parse()`` does the download (if
needed), the read, and the parse in one step, returning the parser so calls can
be chained:

>>> from pidibble.pdbparse import PDBParser
>>> p = PDBParser(source_db='rcsb', source_id='4zmj').parse()

Choosing a source
-----------------

A remote source is selected with two arguments: ``source_db`` names the
database and ``source_id`` is the entry's accession code.  Three databases are
supported:

.. list-table::
   :header-rows: 1
   :widths: 12 20 68

   * - ``source_db``
     - ``source_id`` is…
     - Notes
   * - ``'rcsb'``
     - a 4-character PDB id (e.g. ``'4zmj'``)
     - The RCSB Protein Data Bank.  Fetches ``<id>.pdb`` or, with
       ``input_format='mmCIF'``, ``<id>.cif``.
   * - ``'alphafold'``
     - a UniProt accession (e.g. ``'O46077'``)
     - The AlphaFold model database.  pidibble first fetches the entry metadata
       (JSON), then the predicted-model PDB it points to.
   * - ``'opm'``
     - a PDB id (e.g. ``'7f1r'``)
     - The Orientations of Proteins in Membranes database.  See
       :ref:`opm-dummy-atoms` below.

.. code-block:: python

   # RCSB (the default database for the legacy PDBcode= shorthand)
   p = PDBParser(source_db='rcsb', source_id='4zmj').parse()

   # AlphaFold predicted model for odorant receptor OR2a from D. melanogaster
   p = PDBParser(source_db='alphafold', source_id='O46077').parse()

   # OPM entry (a sweet-taste receptor, with membrane dummy atoms)
   p = PDBParser(source_db='opm', source_id='7f1r').parse()

An unrecognized ``source_db`` raises ``ValueError``, and a ``source_db`` with no
``source_id`` raises ``ValueError`` as well.

The ``PDBcode`` shorthand
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For backward compatibility, ``PDBcode='4zmj'`` is accepted as shorthand for
``source_db='rcsb', source_id='4zmj'``:

>>> p = PDBParser(PDBcode='4zmj').parse()      # doctest: +SKIP

New code should prefer the explicit ``source_db``/``source_id`` form.

Reading a local file
---------------------

To parse a file already on disk, pass ``filepath`` instead of a source:

.. code-block:: python

   p = PDBParser(filepath='mystructure.pdb').parse()
   p = PDBParser(filepath='mystructure.cif', input_format='mmCIF').parse()

If the file does not exist, ``parse()`` raises ``FileNotFoundError``.  No network
access is attempted when ``filepath`` is given.

Caching and re-downloading
--------------------------

For the RCSB and OPM sources, pidibble downloads into the current working
directory and **reuses an existing file** rather than fetching again — so a
second parse of the same entry is offline and instant.  To force a fresh
download (for example, if the local copy is stale or truncated), pass
``overwrite=True``:

.. code-block:: python

   p = PDBParser(source_db='rcsb', source_id='4zmj', overwrite=True).parse()

If a download fails (network error, unknown id), ``fetch()`` logs a warning and
``parse()`` yields an empty :attr:`~pidibble.pdbparse.PDBParser.parsed`.

Choosing the input format
-------------------------

``input_format`` selects the on-disk format and defaults to ``'PDB'``.  Set it to
``'mmCIF'`` to read a PDBx/mmCIF file:

.. code-block:: python

   p = PDBParser(source_db='rcsb', source_id='4tvp', input_format='mmCIF').parse()

The parsed record objects are the same either way; see :doc:`mmcif` for the
coverage and the handful of representational differences.

.. _opm-dummy-atoms:

OPM dummy atoms
---------------

OPM structures embed *dummy* (``DUM``) atoms that mark the position of the
membrane boundary.  When fetching from OPM, pidibble strips blank and ``END``
lines and **splits the dummy atoms into a companion file**: parsing
``7f1r`` writes both ``7f1r.pdb`` (the structure) and ``7f1r-dum.pdb`` (the
membrane markers), so the coordinate parse is not polluted by the dummies.
