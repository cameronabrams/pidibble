.. _quickstart:

Quickstart
==========

This is a five-minute tour.  For the full treatment of each topic, see the
:doc:`User Guide <guide/index>`.

Install pidibble (see :ref:`installation`) and parse an entry.  We'll use
``4ZMJ``, a trimeric ectodomain construct of the HIV-1 envelope glycoprotein:

>>> from pidibble.pdbparse import PDBParser
>>> p = PDBParser(source_db='rcsb', source_id='4zmj').parse()

Constructing the ``PDBParser`` and calling :meth:`~pidibble.pdbparse.PDBParser.parse`
downloads ``4zmj.pdb`` from the RCSB (only if it is not already in the working
directory) and parses it into the dictionary
:attr:`~pidibble.pdbparse.PDBParser.parsed`.

Every parsed record type is a key of ``p.parsed``:

>>> keys = sorted(p.parsed.keys())
>>> len(keys)
60
>>> keys[:8]
['ANISOU', 'ATOM', 'AUTHOR', 'CISPEP', 'COMPND', 'CONECT', 'CRYST1', 'DBREF']

Each value is either a single :class:`~pidibble.pdbrecord.PDBRecord` (for records
that occur once, like ``HEADER``) or a :class:`~pidibble.pdbrecord.PDBRecordList`
of them (for records that recur, like ``ATOM``).  Use
:meth:`~pidibble.baserecord.BaseRecord.pstr` to see what any record holds:

>>> print(p.parsed['HEADER'].pstr())
HEADER
      classification: VIRAL PROTEIN
             depDate: 04-MAY-15
              idCode: 4ZMJ

Fields are plain instance attributes:

>>> p.parsed['HEADER'].classification
'VIRAL PROTEIN'
>>> atoms = p.parsed['ATOM']
>>> len(atoms)
4518
>>> print(atoms[0].pstr())
ATOM
              serial: 1
                name: N
              altLoc:
             residue: resName: LEU; chainID: G; seqNum: 34; iCode:
                   x: -0.092
                   y: 99.33
                   z: 57.967
           occupancy: 1.0
          tempFactor: 137.71
             element: N
              charge:

Notice that structured sub-fields — like the ``residue`` of an atom — are
themselves small record objects with their own attributes:

>>> atoms[0].residue.resName
'LEU'
>>> atoms[0].residue.chainID
'G'
>>> atoms[0].residue.seqNum
34

Reading an mmCIF/PDBx file instead of a legacy PDB file is a one-word change:

>>> c = PDBParser(source_db='rcsb', source_id='4zmj', input_format='mmCIF').parse()  # doctest: +SKIP

pidibble parses the mmCIF entry into the *same* record objects, using author
(``auth``) numbering so the result is a drop-in equivalent of the PDB parse.  See
:doc:`guide/mmcif` for the details and the small, deliberate differences.

Where to go next
----------------

* :doc:`guide/loading` — fetching from the RCSB, AlphaFold and OPM, and reading local files.
* :doc:`guide/data_model` — how ``parsed`` is organized and how to navigate it.
* :doc:`guide/records` — residues, token groups, continuation records, and tables.
* :doc:`guide/assemblies` — biological assemblies and crystallographic symmetry.
* :doc:`guide/mmcif` — PDBx/mmCIF parsing and its coverage.
* :doc:`guide/record_reference` — the full list of supported record types.
