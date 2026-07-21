The parsed data model
=====================

After :meth:`~pidibble.pdbparse.PDBParser.parse`, everything pidibble extracted
lives in one attribute:

>>> from pidibble.pdbparse import PDBParser
>>> p = PDBParser(source_db='rcsb', source_id='4zmj').parse()
>>> p.parsed                                    # doctest: +ELLIPSIS
{...}

``parsed`` is a :class:`~pidibble.pdbrecord.PDBRecordDict` — a dictionary keyed
by record type.

Single records versus record lists
----------------------------------

Each value is one of two things:

* a single :class:`~pidibble.pdbrecord.PDBRecord`, for records that occur once
  per structure (``HEADER``, ``TITLE``, ``CRYST1``, …); or
* a :class:`~pidibble.pdbrecord.PDBRecordList` of ``PDBRecord`` instances, for
  *multiple-entry* records that recur (``ATOM``, ``SEQRES``, ``LINK``, …).

``PDBRecordList`` subclasses :class:`collections.UserList`, so test for it with
``isinstance`` — a plain ``type(v) == list`` check never matches, because a
``UserList`` is not a ``list``:

>>> from pidibble.pdbrecord import PDBRecordList
>>> isinstance(p.parsed['ATOM'], PDBRecordList)      # a multiple-entry record
True
>>> isinstance(p.parsed['HEADER'], PDBRecordList)    # a single-instance record
False
>>> [k for k, v in p.parsed.items() if type(v) == list]   # the naive check never matches
[]

For everyday use, a ``PDBRecordList`` behaves like a list — index it, slice it,
iterate it, take its ``len``:

>>> atoms = p.parsed['ATOM']
>>> len(atoms)
4518
>>> first, last = atoms[0], atoms[-1]

Which records land in which bucket is documented per record type in
:doc:`record_reference`.

Records and their fields
------------------------

A :class:`~pidibble.pdbrecord.PDBRecord` stores each parsed field as an instance
attribute of the appropriate Python type — strings, ``int``, ``float``, or
nested record objects.  The quickest way to see what a record contains is
:meth:`~pidibble.baserecord.BaseRecord.pstr` ("pretty string"):

>>> print(p.parsed['HEADER'].pstr())
HEADER
      classification: VIRAL PROTEIN
             depDate: 04-MAY-15
              idCode: 4ZMJ

The keys shown by ``pstr()`` are exactly the attribute names:

>>> h = p.parsed['HEADER']
>>> h.classification, h.depDate, h.idCode
('VIRAL PROTEIN', '04-MAY-15', '4ZMJ')

``pstr()`` takes an ``excludes`` list (bookkeeping fields ``key``, ``format`` and
``continuation`` are hidden by default) and a ``pad`` width for the labels, so
you can widen or narrow the display or reveal the hidden internals.

Residues and other structured sub-fields
-----------------------------------------

Some fields are not scalars but small structured objects.  The most common is a
**residue**, which groups a residue name, chain id, sequence number and
insertion code:

>>> a = p.parsed['ATOM'][0]
>>> print(a.residue)                            # doctest: +NORMALIZE_WHITESPACE
resName: LEU; chainID: G; seqNum: 34; iCode:
>>> a.residue.resName, a.residue.chainID, a.residue.seqNum
('LEU', 'G', 34)

Records that describe a *relationship between* residues expose more than one:

>>> b = p.parsed['SSBOND'][0]
>>> b.residue1.resName, b.residue1.chainID, b.residue1.seqNum
('CYS', 'G', 54)
>>> b.residue2.resName, b.residue2.chainID, b.residue2.seqNum
('CYS', 'G', 74)

Working with the model as a whole
---------------------------------

Because ``parsed`` is an ordinary mapping, the usual idioms apply — membership
tests, iteration, comprehensions:

>>> 'ATOM' in p.parsed and 'HETATM' in p.parsed
True
>>> n_hetatm = len(p.parsed['HETATM'])
>>> chains = {a.residue.chainID for a in p.parsed['ATOM']}
>>> sorted(chains)
['B', 'G']

A record type that was not present in the file is simply absent from ``parsed``,
so guard with ``in`` (or ``dict.get``) before reaching for an optional record:

>>> 'SPLIT' in p.parsed
False
