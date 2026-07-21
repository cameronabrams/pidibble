Supported record types
=======================

pidibble parses **all standard PDB record types**.  The table below lists them
with the container you get in :attr:`~pidibble.pdbparse.PDBParser.parsed`, whether
the record is also produced from mmCIF input, and a one-line description.  For the
exact fields of any record, see the :doc:`API reference <../api/API>`.

Reading the columns
-------------------

**Container** is what ``p.parsed[<key>]`` holds:

* *single* — one :class:`~pidibble.pdbrecord.PDBRecord`.
* *list* — a :class:`~pidibble.pdbrecord.PDBRecordList` of records (a
  *multiple-entry* record).
* *sub-keys* — the record expands into several keys (e.g. ``JRNL`` →
  ``JRNL.AUTH``, ``JRNL.TITL``, …; ``REMARK`` → ``REMARK.<n>`` and, for tables,
  ``REMARK.350.BIOMOLECULE1.TRANSFORM<n>``).
* *grouping* — a structural marker used while parsing multi-model files.

**mmCIF** marks record types the :doc:`mmcif` path also produces (using author
numbering).  A record type absent from a given file is simply absent from
``parsed``.

.. list-table::
   :header-rows: 1
   :widths: 16 12 8 64

   * - Record
     - Container
     - mmCIF
     - Description
   * - ``ATOM``
     - list
     - ✓
     - Atomic coordinates for standard residues.
   * - ``HETATM``
     - list
     - ✓
     - Atomic coordinates for heterogens (ligands, ions, waters).
   * - ``ANISOU``
     - list
     -
     - Anisotropic temperature factors.
   * - ``TER``
     - list
     -
     - Chain terminator.
   * - ``CONECT``
     - list
     -
     - Explicit atom-to-atom connectivity.
   * - ``LINK``
     - list
     - ✓
     - Covalent and metal-coordination links between residues.
   * - ``SSBOND``
     - list
     - ✓
     - Disulfide bonds.
   * - ``CISPEP``
     - list
     -
     - Cis-peptide bonds.
   * - ``SEQRES``
     - list
     - ✓
     - Primary sequence, one record per chain.
   * - ``DBREF``, ``DBREF1``, ``DBREF2``
     - list
     -
     - Cross-references to sequence databases (``DBREF1``/``DBREF2`` are the
       two-line variant).
   * - ``SEQADV``
     - list
     - ✓
     - Differences between the entry and its sequence database.
   * - ``MODRES``
     - list
     -
     - Modified (non-standard) residues.
   * - ``HELIX``
     - list
     - ✓
     - Helix secondary structure.
   * - ``SHEET``
     - list
     - ✓
     - β-sheet secondary structure (ranges, sense, registration).
   * - ``HET``
     - list
     -
     - Non-standard (heterogen) groups present.
   * - ``HETNAM``
     - list
     -
     - Chemical names of heterogens.
   * - ``HETSYN``
     - list
     -
     - Synonyms for heterogen names.
   * - ``FORMUL``
     - list
     -
     - Chemical formula of each heterogen.
   * - ``SITE``
     - list
     -
     - Functionally important sites.
   * - ``MTRIX1``, ``MTRIX2``, ``MTRIX3``
     - list
     -
     - Non-crystallographic symmetry transforms.
   * - ``REVDAT``
     - list
     -
     - Revision history.
   * - ``HEADER``
     - single
     - ✓
     - Classification, deposition date, and PDB id.
   * - ``TITLE``
     - single
     - ✓
     - Title of the experiment.
   * - ``COMPND``
     - single [#f1]_
     - ✓
     - Description of the macromolecular components (token groups).
   * - ``SOURCE``
     - single [#f1]_
     - ✓
     - Biological source of the components (token groups).
   * - ``KEYWDS``
     - single
     - ✓
     - Keywords describing the entry.
   * - ``EXPDTA``
     - single
     - ✓
     - Experimental technique.
   * - ``AUTHOR``
     - single
     -
     - Depositors.
   * - ``CRYST1``
     - single
     - ✓
     - Unit-cell parameters and space group.
   * - ``ORIGX1``, ``ORIGX2``, ``ORIGX3``
     - single
     -
     - Transformation to the crystallographic origin.
   * - ``SCALE1``, ``SCALE2``, ``SCALE3``
     - single
     -
     - Transformation to fractional crystallographic coordinates.
   * - ``MDLTYP``
     - single
     -
     - Model-type annotations.
   * - ``NUMMDL``
     - single
     -
     - Number of models in the entry.
   * - ``CAVEAT``
     - single
     -
     - Warnings about known errors in the entry.
   * - ``OBSLTE``
     - single
     -
     - Entries made obsolete by this one.
   * - ``SPRSDE``
     - single
     -
     - Entries this one supersedes.
   * - ``SPLIT``
     - single
     -
     - Related entries the structure is split across.
   * - ``MASTER``
     - single
     -
     - Bookkeeping counts of records in the file.
   * - ``END``
     - single
     -
     - End-of-file marker.
   * - ``JRNL``
     - sub-keys
     -
     - Primary citation (``JRNL.AUTH``, ``JRNL.TITL``, ``JRNL.REF``,
       ``JRNL.REFN``, ``JRNL.PMID``, ``JRNL.DOI``).
   * - ``REMARK``
     - sub-keys
     - ✓ [#f2]_
     - Annotation remarks, one key per remark number (``REMARK.<n>``); tables
       expand further.
   * - ``MODEL``, ``ENDMDL``
     - grouping
     -
     - Model boundaries in multi-model entries.

.. rubric:: Notes

.. [#f1] From mmCIF input, ``COMPND`` and ``SOURCE`` are instead flat, per-entity
   :class:`~pidibble.pdbrecord.PDBRecordList`\ s rather than token groups — see
   :doc:`mmcif`.

.. [#f2] For ``REMARK``, the mmCIF path produces the biological-assembly
   transforms (``REMARK 350``) and the missing-residue list (``REMARK 465``).
