PDBx/mmCIF parsing
==================

Many recent RCSB entries have no legacy PDB file at all — the PDBx/mmCIF format
is now the archive's primary format.  pidibble reads mmCIF files and produces the
**same record objects** as the PDB path, so downstream code can stay largely
format-agnostic.

Enabling it
-----------

Pass ``input_format='mmCIF'``.  Everything else — the source arguments, the
:attr:`~pidibble.pdbparse.PDBParser.parsed` dictionary, the record attributes —
is unchanged:

>>> from pidibble.pdbparse import PDBParser
>>> p = PDBParser(source_db='rcsb', source_id='4tvp', input_format='mmCIF').parse()   # doctest: +SKIP

From the RCSB this fetches ``4tvp.cif`` instead of ``4tvp.pdb``.  Under the hood,
pidibble reads the file with the RCSB/wwPDB ``py-mmcif`` library (the
authoritative, C++-backed reader) and maps selected mmCIF *categories* onto PDB
*record types* via a declarative specification.

Author numbering
----------------

pidibble builds the mmCIF result to be the **equivalent of the PDB parse**: it
uses author (``auth_*``) chain ids and sequence numbers, not the internal label
(``label_*``) numbering.  So a residue that a PDB file calls ``chain G, seq 88``
comes out the same way from the mmCIF file, and records like ``LINK`` and
``SSBOND`` line up residue-for-residue with the PDB parse.

What it covers
--------------

As of version 1.7.0 the mmCIF path maps a broad set of record types: coordinates
(``ATOM``/``HETATM``), connectivity (``LINK`` — including metal coordination —
and ``SSBOND``), sequence (``SEQRES``), secondary structure (``HELIX``,
``SHEET``), biological assemblies (``REMARK 350``), missing residues
(``REMARK 465``), sequence-database differences (``SEQADV``), header/metadata
(``HEADER``, ``TITLE``, ``EXPDTA``, ``KEYWDS``, ``CRYST1``) and entity/source
information (``COMPND``/``SOURCE``).  Each mmCIF-derived record mirrors the
attribute names of its PDB counterpart and is validated field-by-field against
the legacy-PDB parse.

The :doc:`record_reference` marks exactly which record types are mmCIF-mapped.
At parse time pidibble also logs, at ``INFO``, how many of the file's mmCIF
categories it read and how many are present but unmapped — a quick way to see
what a given entry carries beyond what pidibble surfaces.

Deliberate departures from PDB equivalence
------------------------------------------

Two differences are intentional and worth knowing about:

``COMPND`` and ``SOURCE`` have a different shape.
   From mmCIF these are emitted as **flat, per-entity records** — a
   :class:`~pidibble.pdbrecord.PDBRecordList` where each element has attributes
   like ``molID``, ``chains`` and ``molName`` (``COMPND``) or ``molID``,
   ``organism``, ``taxid``, ``exprSystem`` and ``gene`` (``SOURCE``) — rather
   than the nested :ref:`token-group <token-groups>` structure the PDB parser
   builds.  The PDB ``FRAGMENT``/``ENGINEERED`` tokens have no clean mmCIF
   equivalent, so **consumers must branch on the input format for these two
   records**.

   .. code-block:: python

      # PDB shape
      p.parsed['COMPND'].tokengroups['compound']['MOL_ID.1'].MOLECULE
      # mmCIF shape
      c.parsed['COMPND'][0].molName

A few purely representational fields differ where the formats themselves differ.
   For example ``HEADER.depDate`` is exposed in mmCIF's native ISO form
   (``2014-06-27``) rather than the PDB ``DD-MON-YY`` form (``27-JUN-14``), and
   ``TITLE``/``KEYWDS`` text is upper-cased to match the PDB convention.

Missing values
--------------

mmCIF writes an inapplicable value as ``.`` and an unknown one as ``?``; an
attribute the file omits is missing too.  All three parse to an empty field, and
because mmCIF carries no per-attribute type, the field's declared type is taken
from the PDB record spec — so a missing value in a field declared ``Float``,
``Integer`` or ``HxInteger`` yields the same guarded
:class:`~pidibble.baseparsers.EmptyField` as it would from a ``.pdb`` file (see
:ref:`empty-fields`).  It equals ``''`` but refuses to act as a number:

.. code-block:: python

   het = cif.parsed['HETATM'][0]
   het.residue.seqNum == ''      # True -- label_seq_id is '.' for a non-polymer
   het.residue.seqNum + 1        # TypeError: Integer field 'seqNum' is empty:
                                 # absent from the mmCIF ('.', '?', or attribute
                                 # not present) ...
   het.residue_auth.seqNum       # 1 -- the auth_* copy carries the number

Comparing the two parses
------------------------

Because the record objects are the same, comparing a PDB and mmCIF parse of the
same entry is straightforward — the important attributes agree, and the mmCIF
values often carry more precision:

.. code-block:: python

   pdb  = PDBParser(source_db='rcsb', source_id='4tvp', input_format='PDB').parse()
   cif  = PDBParser(source_db='rcsb', source_id='4tvp', input_format='mmCIF').parse()

   b_pdb = pdb.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM2']
   b_cif = cif.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM2']
   # same chains (possibly in a different order) and same transform, e.g.
   #   PDB  row t: -515.56   (2 decimals)
   #   mmCIF row t: -515.56  (full precision -0.8660254038 in the rotation)

Current limitations
-------------------

The mmCIF path is deliberately scoped to what downstream consumers need.  It does
not yet group atoms into ``MODEL``/``ENDMDL`` structures for multi-model files,
and it maps a subset of the file's categories (records such as ``ANISOU``,
``CONECT``, ``DBREF``, ``SITE``, ``FORMUL``/``HETNAM`` and ``REVDAT`` are
PDB-only for now).  The :doc:`record_reference` is the authoritative list.

Writing back out
----------------

``write_PDB()`` re-emits ``JRNL`` and ``REMARK`` verbatim from the source
``.pdb`` lines, which an mmCIF parse does not have.  ``JRNL`` is therefore
rebuilt from the ``JRNL.*`` sub-records the ``citation`` category supplies, so
the entry's own paper survives the conversion; see :doc:`citations`.  The
rebuilt block matches the entry's own ``.pdb`` file byte for byte, apart from
one character it cannot know: mmCIF records an ISSN without saying whether it
is the print or the electronic one, so ``REFN`` reads ``ISSN`` where the
``.pdb`` says ``ESSN``.

``REMARK`` is *not* rebuilt.  Most of its content is free text held as an
``LList`` whose line breaks the parser discards, so it cannot be reproduced
faithfully; it is omitted and reported to the logger.

``COMPND`` and ``SOURCE`` are likewise absent from the output.  The mmCIF path
carries them as tokenized attributes (``molID``, ``molName``, ``chains``)
rather than as the raw ``compound``/``srcName`` list the PDB format declares,
so there is nothing in the PDB's own spelling to write.  Rather than emit a
bare card name — which is not a conformant record, and which the tokenizer
cannot re-read — pidibble omits the card and says so at ``INFO``.
