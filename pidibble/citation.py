# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""

.. module:: citation
   :synopsis: normalized bibliographic records for the papers a structure cites

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

A structure file names the paper that describes it, and often several earlier
papers as well. The PDB format files them under ``JRNL`` (the structure's own
paper) and ``REMARK 1 REFERENCE n`` / ``REMARK 0 ORIGINAL DATA REFERENCE``;
mmCIF files all of them in one ``citation`` table.
:func:`~pidibble.citation.citations_from_parsed` and
:func:`~pidibble.citation.citations_from_cif` read either shape and hand back
the same :class:`~pidibble.citation.Citation` objects, so a caller that wants to
credit the work behind a structure never has to know which format it read.

Nothing here touches the network. :func:`~pidibble.citation.citations_from_rcsb`
does, and is reached only when a caller passes ``enrich=True`` to
:meth:`~pidibble.pdbparse.PDBParser.citations`.

"""
from __future__ import annotations

import json
import logging
import re
import urllib.error
import urllib.request
from dataclasses import dataclass, field, asdict

logger = logging.getLogger(__name__)

#: Role a citation plays for the structure that names it.
PRIMARY = 'primary'            #: the paper describing this structure
RELATED = 'related'            #: earlier work the entry cites (REMARK 1 / non-primary ``citation`` rows)
ORIGINAL_DATA = 'original_data'  #: the entry this one re-refines (REMARK 0)

#: ``parsed`` keys that hold citation sub-records, and the role each implies.
_KEY_PATTERNS = (
    (re.compile(r'^JRNL\.(?P<sub>[A-Z]+)$'), PRIMARY, PRIMARY),
    (re.compile(r'^REMARK\.1\.REFERENCE(?P<idx>\d*)\.(?P<sub>[A-Z]+)$'), RELATED, None),
    (re.compile(r'^REMARK\.0\.ORIGREF(?P<idx>\d*)\.(?P<sub>[A-Z]+)$'), ORIGINAL_DATA, None),
)

_ISSN_RE = re.compile(r'^\d{4}-\d{3}[\dxX]$')


def normalize_author(name):
    """
    Render an author name in the ``'Surname, I.N.'`` form every bibliographic
    style wants, from either spelling a structure file uses.

    The PDB format writes ``'M.PANCERA'``; mmCIF writes ``'Pancera, M.'``. The
    latter passes through unchanged.

    Parameters
    ----------
    name : str
        An author name as a structure file spells it.

    Returns
    -------
    str
        The normalized name, e.g. ``'Pancera, M.'``, ``'Bergsten, P.-C.'``.
    """
    name = str(name).strip()
    if not name:
        return ''
    if ',' in name:                 # already mmCIF-style
        last, initials = (x.strip() for x in name.split(',', 1))
        return f'{last}, {initials}' if initials else last
    cut = name.rfind('.')
    if cut < 0:                     # a bare surname, no initials
        return name.title() if name.isupper() else name
    initials, last = name[:cut + 1], name[cut + 1:]
    last = last.title() if last.isupper() else last
    return f'{last}, {initials}' if last else initials


def _as_record(value):
    """``parsed`` holds either a record or a list of them; take the first."""
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value


def _get(rec, attr, default=''):
    """Read one field, treating a guarded empty (``EmptyField``) as absent."""
    if rec is None:
        return default
    val = getattr(rec, attr, default)
    if val is None or val == '':
        return default
    # EmptyField stands in for a numeric field the source left blank
    if type(val).__name__ == 'EmptyField':
        return default
    return val


def _as_int(val):
    """Coerce to int, or None if the value is absent or not a number."""
    try:
        return int(str(val).strip())
    except (TypeError, ValueError):
        return None


def _as_pmid(val):
    """Coerce to a PubMed ID. mmCIF and the RCSB API spell 'this paper has no
    PubMed ID' as ``-1`` (1MBN, for one), which is not an identifier."""
    pmid = _as_int(val)
    return pmid if (pmid is not None and pmid > 0) else None


def _clean_orcids(values):
    """Keep only real ORCIDs -- both sources pad the list with '?' placeholders
    so it stays parallel to the author list."""
    return [str(v) for v in (values or []) if str(v).strip() not in ('', '.', '?')]


@dataclass
class Citation:
    """
    One paper a structure entry credits, in a form independent of the file
    format it was read from.

    ``authors`` is always ``'Surname, I.N.'``; ``authors_raw`` preserves the
    file's own spelling. ``title`` is verbatim -- mmCIF carries proper case,
    the PDB format carries upper case, and pidibble does not guess sentence
    case back because doing so mangles acronyms.

    Case is the one thing that does not survive a PDB-format source: it stores
    ``'A.B.MCDERMOTT'``, from which ``'McDermott, A.B.'`` is not recoverable
    (nor ``'de Val'`` from ``'N.DE VAL'``), so those normalize to
    ``'Mcdermott, A.B.'`` and ``'De Val, N.'``. Read the mmCIF, or pass
    ``enrich=True``, wherever the exact rendering matters -- both carry the
    author's own capitalization. Identifiers, years, volumes and pages are
    unaffected and agree exactly across the two formats.
    """
    pdb_id: str = ''
    id: str = ''
    role: str = PRIMARY
    title: str = ''
    authors: list = field(default_factory=list)
    authors_raw: list = field(default_factory=list)
    editors: list = field(default_factory=list)
    orcids: list = field(default_factory=list)
    journal: str = ''
    journal_full: str = ''
    publisher: str = ''
    volume: str = ''
    first_page: str = ''
    last_page: str = ''
    year: int | None = None
    doi: str | None = None
    pmid: int | None = None
    issn: str | None = None
    isbn: str | None = None
    source: str = ''
    raw: dict = field(default_factory=dict)

    @property
    def doi_url(self):
        """Resolvable ``https://doi.org/...`` URL, or None if there is no DOI."""
        return f'https://doi.org/{self.doi}' if self.doi else None

    @property
    def pubmed_url(self):
        """PubMed URL for this paper, or None if there is no PMID."""
        return f'https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/' if self.pmid else None

    @property
    def is_primary(self):
        """True for the paper describing this structure."""
        return self.role == PRIMARY

    def asdict(self):
        """Return this citation as a plain dict, for callers that want no
        dependency on the dataclass."""
        return asdict(self)

    def reference_string(self):
        """
        A single-line human-readable reference, for a "please cite" report.

        Returns
        -------
        str
            e.g. ``'Pancera, M., Zhou, T., ... Kwong, P.D. (2014) Structure and
            immune recognition of trimeric pre-fusion HIV-1 Env. Nature
            514:455-461. doi:10.1038/nature13808'``
        """
        parts = []
        if self.authors:
            parts.append(', '.join(self.authors))
        if self.year:
            parts.append(f'({self.year})')
        if self.title:
            parts.append(self.title if self.title.endswith('.') else self.title + '.')
        venue = self.journal or self.journal_full
        if venue:
            pages = self.first_page
            if self.last_page and self.last_page != self.first_page:
                pages = f'{self.first_page}-{self.last_page}'
            loc = f'{venue} {self.volume}:{pages}' if (self.volume and pages) else venue
            parts.append(loc.rstrip(':') + '.')
        if self.publisher:
            parts.append(self.publisher + '.')
        if self.doi:
            parts.append(f'doi:{self.doi}')
        return ' '.join(p for p in parts if p)

    def __str__(self):
        return self.reference_string()


@dataclass
class CitationId:
    """
    The bare identifiers for a paper -- enough to name it unambiguously and
    hand it to a resolver, with none of the reconstruction that a formatted
    reference string involves.
    """
    pdb_id: str = ''
    id: str = ''
    role: str = PRIMARY
    doi: str | None = None
    pmid: int | None = None

    @property
    def doi_url(self):
        """Resolvable ``https://doi.org/...`` URL, or None if there is no DOI."""
        return f'https://doi.org/{self.doi}' if self.doi else None

    @property
    def pubmed_url(self):
        """PubMed URL for this paper, or None if there is no PMID."""
        return f'https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/' if self.pmid else None

    def asdict(self):
        """Return these identifiers as a plain dict."""
        return asdict(self)


def _group_citation_keys(parsed):
    """
    Bucket the citation sub-records in ``parsed`` by the citation they belong to.

    Returns
    -------
    dict
        ``{(role, id): {subkey: record}}``, in the order the keys were found.
    """
    groups = {}
    for key, value in parsed.items():
        for pattern, role, fixed_id in _KEY_PATTERNS:
            m = pattern.match(key)
            if not m:
                continue
            cid = fixed_id if fixed_id is not None else (m.groupdict().get('idx') or '1')
            groups.setdefault((role, cid), {})[m.group('sub')] = _as_record(value)
            break
    return groups


def _citation_from_subrecords(subs, pdb_id, role, cid):
    """Assemble one :class:`Citation` from a ``{subkey: record}`` bucket."""
    ref = subs.get('REF')
    refn = subs.get('REFN')
    raw_authors = list(_get(subs.get('AUTH'), 'authorList', []) or [])
    raw_editors = list(_get(subs.get('EDIT'), 'editorList', []) or [])
    serial = str(_get(refn, 'issn', ''))
    kind = str(_get(refn, 'issnORessn', '')).upper()
    doi = str(_get(subs.get('DOI'), 'doi', '')).lower() or None
    return Citation(
        pdb_id=pdb_id,
        id=cid,
        role=role,
        title=str(_get(subs.get('TITL'), 'title', '')),
        authors=[normalize_author(a) for a in raw_authors],
        authors_raw=raw_authors,
        editors=[normalize_author(e) for e in raw_editors],
        journal=str(_get(ref, 'pubName', '')),
        publisher=str(_get(subs.get('PUBL'), 'pub', '')),
        volume=str(_get(ref, 'volume', '')),
        first_page=str(_get(ref, 'page', '')),
        year=_as_int(_get(ref, 'year', None)),
        doi=doi,
        pmid=_as_pmid(_get(subs.get('PMID'), 'pmid', None)),
        issn=serial if (kind != 'ISBN' and serial) else None,
        isbn=serial if (kind == 'ISBN' and serial) else None,
        source='JRNL' if role == PRIMARY else 'REMARK1',
        raw={'subrecords': sorted(subs), 'doi_raw': str(_get(subs.get('DOI'), 'doi', '')) or None},
    )


def citations_from_parsed(parsed, pdb_id=''):
    """
    Build citations from parsed records, whichever format produced them.

    This reads only ``JRNL.*``, ``REMARK.1.REFERENCE<n>.*`` and
    ``REMARK.0.ORIGREF<n>.*`` -- the keys an mmCIF parse emits too -- so it
    works on both input paths. :func:`citations_from_cif` gets more out of an
    mmCIF file (proper title case, last page, verbatim author spelling) and is
    preferred where the ``citation`` category is at hand.

    Parameters
    ----------
    parsed : PDBRecordDict
        A parser's :attr:`~pidibble.pdbparse.PDBParser.parsed` dictionary.
    pdb_id : str, optional
        Stamped onto every citation so a caller can tell whose paper it is.

    Returns
    -------
    list of Citation
        Primary first, then related, then original-data references.
    """
    cits = [_citation_from_subrecords(subs, pdb_id, role, cid)
            for (role, cid), subs in _group_citation_keys(parsed).items()]
    return _ordered(cits)


def citations_from_cif(cif_data, pdb_id=''):
    """
    Build citations straight from the mmCIF ``citation``/``citation_author``
    categories, keeping everything the PDB-format layout cannot hold.

    Parameters
    ----------
    cif_data : DataContainer
        The container a parser holds in :attr:`~pidibble.pdbparse.PDBParser.cif_data`.
    pdb_id : str, optional
        Stamped onto every citation.

    Returns
    -------
    list of Citation
        Primary first, then the rest in file order. Empty if the file carries
        no ``citation`` category.
    """
    if cif_data is None:
        return []
    cit = cif_data.getObj('citation')
    if cit is None or not len(cit):
        return []

    auth = cif_data.getObj('citation_author')
    by_id = {}
    orcid_by_id = {}
    if auth is not None:
        rows = [auth.getRowAttributeDict(i) for i in range(len(auth))]

        def ordinal(row):
            try:
                return int(row.get('ordinal', 0))
            except (TypeError, ValueError):
                return 0

        for row in sorted(rows, key=ordinal):
            name = str(row.get('name', '')).strip()
            if name and name not in '.?':
                cid = str(row.get('citation_id', '')).strip()
                by_id.setdefault(cid, []).append(name)
                orcid_by_id.setdefault(cid, []).append(row.get('identifier_ORCID', ''))

    def val(row, attr):
        v = row.get(attr, '')
        v = '' if v is None else str(v).strip()
        return '' if v in ('.', '?') else v

    out = []
    nonprimary = 0
    for i in range(len(cit)):
        row = cit.getRowAttributeDict(i)
        cid = str(row.get('id', '')).strip()
        if cid == 'primary':
            role, label = PRIMARY, PRIMARY
        else:
            nonprimary += 1
            role, label = RELATED, (cid if cid.isdigit() else str(nonprimary))
        raw_authors = by_id.get(cid, [])
        isbn, issn = val(row, 'book_id_ISBN'), val(row, 'journal_id_ISSN')
        serial = isbn or issn
        is_isbn = bool(serial) and (bool(isbn) or not _ISSN_RE.match(serial))
        doi = val(row, 'pdbx_database_id_DOI').lower() or None
        out.append(Citation(
            pdb_id=pdb_id,
            id=label,
            role=role,
            title=val(row, 'title'),
            authors=[normalize_author(a) for a in raw_authors],
            authors_raw=list(raw_authors),
            orcids=_clean_orcids(orcid_by_id.get(cid)),
            journal=val(row, 'journal_abbrev') or val(row, 'book_title'),
            journal_full=val(row, 'journal_full'),
            publisher=val(row, 'book_publisher'),
            volume=val(row, 'journal_volume'),
            first_page=val(row, 'page_first'),
            last_page=val(row, 'page_last'),
            year=_as_int(val(row, 'year')),
            doi=doi,
            pmid=_as_pmid(val(row, 'pdbx_database_id_PubMed')),
            issn=None if is_isbn else (serial or None),
            isbn=serial if is_isbn else None,
            source='mmCIF',
            raw={'citation_id': cid, 'doi_raw': val(row, 'pdbx_database_id_DOI') or None},
        ))
    return _ordered(out)


#: Where :func:`citations_from_rcsb` looks. The entry endpoint returns the same
#: citation data the RCSB entry page shows.
RCSB_ENTRY_URL = 'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}'


def citations_from_rcsb(pdb_id, timeout=10, url_template=RCSB_ENTRY_URL):
    """
    Fetch citations for a PDB ID from the RCSB Data API -- the same records the
    RCSB entry page displays.

    **This contacts the network.** Nothing in pidibble calls it unless a caller
    asks for it explicitly via ``citations(enrich=True)``.

    The API is worth the round trip when the local file is a ``.pdb``: it
    returns properly-cased titles and journal names, the last page, and author
    names already in ``'Surname, I.N.'`` form -- none of which the PDB format
    carries.

    Parameters
    ----------
    pdb_id : str
        A four-character PDB ID.
    timeout : float, optional
        Socket timeout in seconds (default 10).
    url_template : str, optional
        Overridable for testing; must contain ``{pdb_id}``.

    Returns
    -------
    list of Citation
        Empty if the entry has no citations or the request failed -- a failure
        is logged, never raised, so an enrichment attempt cannot break a parse.
    """
    if not pdb_id:
        return []
    url = url_template.format(pdb_id=pdb_id.lower())
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            payload = json.load(response)
    except (urllib.error.URLError, OSError, ValueError) as e:
        logger.warning(f'Could not fetch citations for {pdb_id} from the RCSB Data API: {e}')
        return []

    # `citation` carries every paper and flags the primary one; the separate
    # `rcsb_primary_citation` object repeats it, and is the only source when an
    # entry has no `citation` list at all.
    entries = list(payload.get('citation') or [])
    if not entries and payload.get('rcsb_primary_citation'):
        entries = [dict(payload['rcsb_primary_citation'], rcsb_is_primary='Y')]

    out = []
    nonprimary = 0
    for entry in entries:
        cid = str(entry.get('id', ''))
        if entry.get('rcsb_is_primary') == 'Y' or cid == PRIMARY:
            role, label = PRIMARY, PRIMARY
        else:
            nonprimary += 1
            role, label = RELATED, (cid if cid.isdigit() else str(nonprimary))
        raw_authors = [str(a) for a in (entry.get('rcsb_authors') or [])]
        doi = str(entry.get('pdbx_database_id_DOI', '') or '').lower() or None
        isbn = str(entry.get('book_id_ISBN', '') or '')
        issn = str(entry.get('journal_id_ISSN', '') or '')
        serial = isbn or issn
        is_isbn = bool(serial) and (bool(isbn) or not _ISSN_RE.match(serial))
        out.append(Citation(
            pdb_id=pdb_id,
            id=label,
            role=role,
            title=str(entry.get('title', '') or ''),
            authors=[normalize_author(a) for a in raw_authors],
            authors_raw=raw_authors,
            orcids=_clean_orcids(entry.get('rcsb_ORCID_identifiers')),
            journal=str(entry.get('rcsb_journal_abbrev') or entry.get('journal_abbrev', '') or ''),
            journal_full=str(entry.get('journal_full', '') or ''),
            publisher=str(entry.get('book_publisher', '') or ''),
            volume=str(entry.get('journal_volume', '') or ''),
            first_page=str(entry.get('page_first', '') or ''),
            last_page=str(entry.get('page_last', '') or ''),
            year=_as_int(entry.get('year')),
            doi=doi,
            pmid=_as_pmid(entry.get('pdbx_database_id_PubMed')),
            issn=None if is_isbn else (serial or None),
            isbn=serial if is_isbn else None,
            source='rcsb-api',
            raw={'citation_id': cid},
        ))
    return _ordered(out)


_ROLE_ORDER = {PRIMARY: 0, RELATED: 1, ORIGINAL_DATA: 2}


def _ordered(citations):
    """Primary first, then related, then original-data, numbered citations in
    numeric order within a role.

    mmCIF does not oblige its ``citation`` rows to appear in id order -- 4HHB
    files them 1, 3, 4, 5, 6, 2, 7, 8 -- while REMARK 1 always counts up, so
    sorting is what makes the two paths agree.
    """
    def key(c):
        cid = str(c.id)
        return (_ROLE_ORDER.get(c.role, 9), int(cid) if cid.isdigit() else 1 << 30, cid)

    return sorted(citations, key=key)


def _merge(base, better):
    """
    Overlay ``better`` onto ``base`` field by field, keeping whatever ``base``
    has where ``better`` is silent.

    Used to let an RCSB Data API answer fill in what the PDB format cannot
    carry without discarding anything already read from the file.
    """
    merged = Citation(**asdict(base))
    for name, value in asdict(better).items():
        if name in ('raw', 'source', 'id', 'role', 'pdb_id'):
            continue
        if value not in ('', None, []):
            setattr(merged, name, value)
    merged.source = f'{base.source}+{better.source}'
    merged.raw = {**base.raw, **{f'rcsb_{k}': v for k, v in better.raw.items()}}
    return merged


def merge_citations(local, remote):
    """
    Combine locally-parsed citations with a richer remote answer.

    Matching is by ``(role, id)``; a remote citation with no local counterpart
    is appended.

    Parameters
    ----------
    local : list of Citation
        Citations read from the file.
    remote : list of Citation
        Citations from :func:`citations_from_rcsb`.

    Returns
    -------
    list of Citation
    """
    if not remote:
        return list(local)
    by_key = {(r.role, r.id): r for r in remote}
    out = []
    for c in local:
        match = by_key.get((c.role, c.id))
        out.append(_merge(c, match) if match else c)
    seen = {(c.role, c.id) for c in local}
    out.extend(r for r in remote if (r.role, r.id) not in seen)
    return _ordered(out)
