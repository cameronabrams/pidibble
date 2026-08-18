Citations: the papers behind a structure
========================================

Every deposited structure names the paper that describes it, and many name
earlier work as well.  Pidibble reads those and hands them back as normalized
:class:`~pidibble.citation.Citation` objects, so a tool that wants to credit the
science behind a model — a "please cite" report, a bibliography, a reference
manager — never has to know which file format it read.

.. code-block:: python

   from pidibble.pdbparse import PDBParser

   p = PDBParser(source_db='rcsb', source_id='4tvp').parse()
   for c in p.citations():
       print(c.reference_string())

.. code-block:: text

   Pancera, M., Zhou, T., Druz, A., ... Kwong, P.D. (2014) STRUCTURE AND IMMUNE
   RECOGNITION OF TRIMERIC PRE-FUSION HIV-1 ENV. NATURE 514:455. doi:10.1038/nature13808

The shouting is the PDB format's, not pidibble's — it stores titles and journal
names in upper case, and pidibble hands them back as it finds them.  Read the
same entry as mmCIF and you get ``Structure and immune recognition of trimeric
pre-fusion HIV-1 Env. Nature 514:455-461``.  See `Author names`_ below.

Where the data comes from
-------------------------

The two file formats file the same information differently, and pidibble reads
both into the same shape:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Role
     - PDB format
     - mmCIF
   * - ``'primary'``
     - ``JRNL``
     - ``citation`` row with ``id == 'primary'``
   * - ``'related'``
     - ``REMARK 1 REFERENCE n``
     - every other ``citation`` row
   * - ``'original_data'``
     - ``REMARK 0 ORIGINAL DATA REFERENCE``
     - —

An mmCIF parse also *emits* the ``JRNL.*`` and ``REMARK.1.REFERENCE<n>.*``
records, so code that reads :attr:`~pidibble.pdbparse.PDBParser.parsed` directly
works unchanged on either input.

Roles
-----

``citations()`` returns everything, primary first.  Filter with ``role``:

.. code-block:: python

   paper   = p.citations(role='primary')[0]      # what a user of this structure owes
   earlier = p.citations(role='related')         # prior work the entry cites

For a citation report, the primary paper is normally what you want; related
references are the kind of thing to put behind an opt-in.

Just the identifiers
--------------------

A DOI and a PMID name a paper unambiguously, and cannot be wrong the way a
reconstructed reference string can.  If your code resolves the bibliographic
detail itself — from Crossref, PubMed, or an existing reference database — ask
for the identifiers alone:

.. code-block:: python

   >>> p.citation_ids()
   [CitationId(pdb_id='4tvp', id='primary', role='primary',
               doi='10.1038/nature13808', pmid=25296255)]

Citations carrying neither identifier are omitted — there is nothing to resolve
them by.  :attr:`~pidibble.citation.CitationId.doi_url` and
:attr:`~pidibble.citation.CitationId.pubmed_url` give resolvable links.

Author names
------------

The PDB format writes ``'M.PANCERA'``; mmCIF writes ``'Pancera, M.'``.
:attr:`~pidibble.citation.Citation.authors` is always the second form — what
BibTeX and every downstream style want — and
:attr:`~pidibble.citation.Citation.authors_raw` preserves the file's own
spelling.

Nothing else about the record depends on which format you read.  Capitalization
does: a PDB-format file stores ``'A.B.MCDERMOTT'``, and ``'McDermott, A.B.'`` is
not recoverable from that, so it normalizes to ``'Mcdermott, A.B.'``.  Likewise
titles, which the PDB format stores in upper case — pidibble hands them back
verbatim rather than guessing sentence case, because guessing mangles acronyms
(``HIV-1``, ``A``, ``SCN``).  Read the mmCIF, or enrich (below), wherever the
exact rendering matters.  Identifiers, years, volumes and pages agree exactly
either way.

Nothing here touches the network
--------------------------------

:meth:`~pidibble.pdbparse.PDBParser.citations` and
:meth:`~pidibble.pdbparse.PDBParser.citation_ids` read the file you already
have.  That is deliberate: a build that silently contacted a web service to
assemble a citation list would be a surprising thing for a citation list to do.

When you *do* want more than the file carries, ask for it:

.. code-block:: python

   rich = p.citations(enrich=True)

This queries the `RCSB Data API <https://data.rcsb.org/>`_ — the same records
the RCSB entry page shows — and overlays properly-cased titles and journal
names, the last page, ORCIDs, and author names in the author's own
capitalization, without discarding anything read from the file.  A failed
request is logged and ignored: enrichment can add to a citation list, never
break one.

Enrichment needs a PDB ID.  It uses ``source_id`` if the parser has one, and
otherwise falls back to the ``HEADER`` record's ``idCode``.

Building your own reference strings
-----------------------------------

:meth:`~pidibble.citation.Citation.reference_string` produces a single-line
reference suitable for a console report.  For anything more structured, the
fields are all there:

.. code-block:: python

   c = p.citations(role='primary')[0]
   c.title, c.authors, c.journal, c.volume, c.first_page, c.last_page, c.year
   c.doi, c.pmid, c.issn, c.isbn, c.publisher, c.editors, c.orcids
   c.doi_url, c.pubmed_url

:meth:`~pidibble.citation.Citation.asdict` returns the whole thing as a plain
dict for callers that would rather not depend on the dataclass.
