Working with records
=====================

:doc:`data_model` covered scalar fields and residue sub-fields.  A few record
types carry richer structure — token groups, grouped lists, continuation text
and tables.  This page shows how to get at each.

>>> from pidibble.pdbparse import PDBParser
>>> p = PDBParser(source_db='rcsb', source_id='4zmj').parse()

.. _token-groups:

Token groups: ``COMPND`` and ``SOURCE``
----------------------------------------

``COMPND`` and ``SOURCE`` are *token-group* records: a single record whose text
is parsed into groups of key/value tokens, one group per molecule (``MOL_ID``).
The parsed groups live in the ``tokengroups`` attribute, keyed first by field
name and then by ``MOL_ID.N``:

>>> compnd = p.parsed['COMPND']
>>> len(compnd.compound)                 # tokens across all molecules
10
>>> mol1 = compnd.tokengroups['compound']['MOL_ID.1']
>>> mol1.MOLECULE
'ENVELOPE GLYCOPROTEIN GP160'
>>> mol1.CHAIN
['G']
>>> mol1.MUTATION
'YES'

``SOURCE`` works the same way (its field is called ``srcName``):

>>> src1 = p.parsed['SOURCE'].tokengroups['srcName']['MOL_ID.1']
>>> src1.ORGANISM_SCIENTIFIC
'HUMAN IMMUNODEFICIENCY VIRUS 1'

Rather than walk the groups by hand, :meth:`~pidibble.pdbrecord.PDBRecord.get_token`
collects one token across every group.  It returns a single value if only one
group has the token, or a ``{MOL_ID.N: value}`` dict if several do:

>>> p.parsed['COMPND'].get_token('CHAIN')
{'MOL_ID.1': ['G'], 'MOL_ID.2': ['B']}

.. note::

   When ``COMPND`` and ``SOURCE`` come from an mmCIF file they are shaped
   differently — flat per-entity records rather than token groups — because the
   PDB tokens have no clean mmCIF equivalent.  See :doc:`mmcif` for that shape and
   why consumers must branch on the input format for these two records.

Grouped records: ``SEQRES``
---------------------------

``SEQRES`` is a *multiple-entry* record: one :class:`~pidibble.pdbrecord.PDBRecord`
per chain, each already assembled from the file's several continuation lines into
a single residue-name list:

>>> seqres = p.parsed['SEQRES']              # a PDBRecordList, one per chain
>>> chain_g = seqres[0]
>>> chain_g.chainID, chain_g.numRes
('G', 481)
>>> chain_g.resNames[:5]
['ALA', 'GLU', 'ASN', 'LEU', 'TRP']
>>> len(chain_g.resNames)
481

Records that relate residues
-----------------------------

Many record types describe a relationship and expose the participating residues
as named sub-fields.  ``LINK`` (covalent/coordination bonds) names two atoms and
their residues:

>>> link = p.parsed['LINK'][0]
>>> print(link.pstr())
LINK
               name1: ND2
             altLoc1:
            residue1: resName: ASN; chainID: G; seqNum: 156; iCode:
               name2: C1
             altLoc2:
            residue2: resName: NAG; chainID: G; seqNum: 615; iCode:
                sym1: 1555
                sym2: 1555
              length: 1.44

``CISPEP`` (cis-peptide bonds) is similar, with numeric fields parsed to their
native types:

>>> cis = p.parsed['CISPEP'][0]
>>> cis.residue1.resName, cis.residue2.resName
('ILE', 'THR')
>>> cis.measure                              # the omega angle, a float
-24.5

Continuation records
--------------------

Several PDB record types span multiple physical lines using a continuation
column.  pidibble stitches these back together for you, so you never see the raw
continuation numbering: multi-line free text (as in ``TITLE``) arrives as one
joined string, and repeated data (as in ``SEQRES`` above) arrives as one list.
The bookkeeping fields that drove the stitching (``continuation`` and friends)
are hidden by :meth:`~pidibble.baserecord.BaseRecord.pstr` by default.

Tables
------

A handful of records embed *tables* — most importantly the ``BIOMT``
transformation matrices under ``REMARK 350`` and the crystallographic symmetry
operators under ``REMARK 290``.  pidibble parses these into their own indexed
sub-records (``REMARK.350.BIOMOLECULE1.TRANSFORM1`` and so on).  Because these
tables are how biological assemblies and symmetry are represented, they get
their own page: :doc:`assemblies`.
