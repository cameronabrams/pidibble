"""

.. module:: test_citation
   :synopsis: the papers a structure entry credits, read from either file format

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

Everything here runs offline against the fixtures in ``test_rcsb/``. The one
test that exercises the RCSB Data API stubs the HTTP call.

"""
import io
import json
import urllib.error
import urllib.request

import pytest

from pidibble.citation import (Citation, CitationId, PRIMARY, RELATED,
                               citations_from_rcsb, merge_citations,
                               normalize_author)
from pidibble.pdbparse import PDBParser

PDB_1CA2 = 'test_rcsb/1ca2.pdb'
CIF_1CA2 = 'test_rcsb/1ca2.cif'
PDB_4TVP = 'test_rcsb/4tvp.pdb'
CIF_4TVP = 'test_rcsb/4tvp.cif'


def parse(path, fmt='PDB'):
    return PDBParser(input_format=fmt, filepath=path).parse()


# --- author normalization --------------------------------------------------

@pytest.mark.parametrize('given,expected', [
    ('M.PANCERA', 'Pancera, M.'),
    ('A.E.ERIKSSON', 'Eriksson, A.E.'),
    ('P.-C.BERGSTEN', 'Bergsten, P.-C.'),
    ('G.B.STEWART-JONES', 'Stewart-Jones, G.B.'),
    ('Pancera, M.', 'Pancera, M.'),          # mmCIF form passes through
    ('Bergsten, P.-C.', 'Bergsten, P.-C.'),
    ('LILJAS', 'Liljas'),                    # surname with no initials
    ('', ''),
])
def test_author_names_normalize_to_one_form(given, expected):
    assert normalize_author(given) == expected


# --- the two input paths converge -----------------------------------------

def test_pdb_and_mmcif_yield_the_same_citations():
    """The same entry read either way credits the same papers, with the same
    authors, identifiers and roles."""
    from_pdb = parse(PDB_1CA2).citations()
    from_cif = parse(CIF_1CA2, 'mmCIF').citations()
    assert len(from_pdb) == len(from_cif) == 5   # primary + REMARK 1 REFERENCE 1-4
    for a, b in zip(from_pdb, from_cif):
        assert (a.role, a.id) == (b.role, b.id)
        assert a.authors == b.authors
        assert a.doi == b.doi
        assert a.pmid == b.pmid
        assert a.year == b.year
        assert a.volume == b.volume
        assert a.first_page == b.first_page


def test_primary_citation_fields():
    c = parse(PDB_1CA2).citations(role=PRIMARY)[0]
    assert c.is_primary
    assert c.pdb_id == '1CA2'                     # falls back to the HEADER idCode
    assert c.authors == ['Eriksson, A.E.', 'Jones, T.A.', 'Liljas, A.']
    assert c.authors_raw == ['A.E.ERIKSSON', 'T.A.JONES', 'A.LILJAS']
    assert c.journal == 'PROTEINS'
    assert (c.volume, c.first_page, c.year) == ('4', '274', 1988)
    assert c.doi == '10.1002/prot.340040406'      # canonical lower case
    assert c.raw['doi_raw'] == '10.1002/PROT.340040406'
    assert c.pmid == 3151019
    assert c.issn == '0887-3585'
    assert c.doi_url == 'https://doi.org/10.1002/prot.340040406'
    assert c.pubmed_url == 'https://pubmed.ncbi.nlm.nih.gov/3151019/'
    assert c.source == 'JRNL'


def test_mmcif_carries_what_the_pdb_format_cannot():
    c = parse(CIF_1CA2, 'mmCIF').citations(role=PRIMARY)[0]
    assert c.title == 'Refined structure of human carbonic anhydrase II at 2.0 A resolution.'
    assert c.last_page == '282'                   # no slot for this in JRNL REF
    assert c.source == 'mmCIF'


def test_pdb_titles_stay_upper_case():
    """pidibble does not guess sentence case back from the PDB format's upper
    case -- that mangles acronyms. It hands back what the file says."""
    c = parse(PDB_1CA2).citations(role=PRIMARY)[0]
    assert c.title == 'REFINED STRUCTURE OF HUMAN CARBONIC ANHYDRASE II AT 2.0 A RESOLUTION.'


# --- REMARK 1 holds more than one reference --------------------------------

def test_each_remark_1_reference_is_its_own_citation():
    """REMARK 1 runs REFERENCE 1, 2, 3 ... back to back with no blank line
    between them; each is a separate paper, not one merged blob."""
    related = parse(PDB_1CA2).citations(role=RELATED)
    assert [c.id for c in related] == ['1', '2', '3', '4']
    assert [c.year for c in related] == [1988, 1975, 1972, 1972]
    assert related[0].authors == ['Eriksson, A.E.', 'Kylsten, P.M.', 'Jones, T.A.', 'Liljas, A.']
    assert related[1].authors == ['Notstrand, B.', 'Vaara, I.', 'Kannan, K.K.']
    assert related[3].journal == 'NATURE NEW BIOL.'
    # each title belongs to exactly one reference
    assert related[3].title == 'CRYSTAL STRUCTURE OF HUMAN CARBONIC ANHYDRASE C'
    assert 'CARBONIC ANHYDRASE C' not in related[0].title


def test_remark_1_reference_records_are_keyed_separately():
    parsed = parse(PDB_1CA2).parsed
    for n in (1, 2, 3, 4):
        assert f'REMARK.1.REFERENCE{n}.TITL' in parsed
    assert parsed['REMARK.1.REFERENCE2.TITL'].title.startswith('STRUCTURAL RELATIONSHIP')


def test_book_reference_carries_editor_and_publisher():
    c = parse(PDB_1CA2).citations(role=RELATED)[1]
    assert c.editors == ['Markert, C.L.']
    assert c.publisher == 'ACADEMIC PRESS,NEW YORK'


# --- mmCIF-sourced records mirror the PDB layout ---------------------------

def test_mmcif_emits_the_pdb_citation_record_keys():
    """A CIF parse must expose JRNL.* and REMARK.1.REFERENCE<n>.* too, so code
    reading records rather than citations works on both paths."""
    parsed = parse(CIF_1CA2, 'mmCIF').parsed
    assert parsed['JRNL.AUTH'].authorList == ['A.E.ERIKSSON', 'T.A.JONES', 'A.LILJAS']
    assert parsed['JRNL.DOI'].doi == '10.1002/PROT.340040406'
    assert parsed['JRNL.PMID'].pmid == 3151019
    assert parsed['JRNL.REF'].year == 1988
    for n in (1, 2, 3, 4):
        assert f'REMARK.1.REFERENCE{n}.AUTH' in parsed


def test_mmcif_citation_records_match_the_pdb_ones():
    from_pdb = parse(PDB_1CA2).parsed
    from_cif = parse(CIF_1CA2, 'mmCIF').parsed
    for key in ('JRNL.AUTH', 'JRNL.PMID', 'JRNL.DOI', 'JRNL.TITL'):
        for attr in ('authorList', 'pmid', 'doi', 'title'):
            if hasattr(from_pdb[key], attr):
                assert getattr(from_cif[key], attr) == getattr(from_pdb[key], attr), key


def test_entry_with_a_single_citation():
    for path, fmt in ((PDB_4TVP, 'PDB'), (CIF_4TVP, 'mmCIF')):
        cits = parse(path, fmt).citations()
        assert len(cits) == 1
        assert cits[0].role == PRIMARY
        assert cits[0].pmid == 25296255
        assert cits[0].doi == '10.1038/nature13808'
        assert len(cits[0].authors) == 27
        assert cits[0].authors[0] == 'Pancera, M.'
        assert cits[0].authors[-1] == 'Kwong, P.D.'


# --- the identifier fast path ----------------------------------------------

def test_citation_ids_returns_only_resolvable_papers():
    ids = parse(PDB_1CA2).citation_ids()
    assert [type(i) for i in ids] == [CitationId]
    # only the primary paper of 1CA2 carries a DOI or PMID
    assert len(ids) == 1
    assert ids[0].role == PRIMARY
    assert ids[0].doi == '10.1002/prot.340040406'
    assert ids[0].pmid == 3151019
    assert ids[0].pubmed_url.endswith('/3151019/')


def test_citation_ids_agree_across_formats():
    assert ([i.asdict() for i in parse(PDB_1CA2).citation_ids()]
            == [i.asdict() for i in parse(CIF_1CA2, 'mmCIF').citation_ids()])


def test_role_filter():
    p = parse(PDB_1CA2)
    assert all(c.role == RELATED for c in p.citations(role=RELATED))
    assert p.citations(role='original_data') == []


# --- reference strings -----------------------------------------------------

def test_reference_string_reads_as_a_reference():
    c = parse(CIF_4TVP, 'mmCIF').citations()[0]
    s = c.reference_string()
    assert s.startswith('Pancera, M., Zhou, T., Druz, A.')
    assert '(2014)' in s
    assert 'Nature 514:455-461.' in s
    assert s.endswith('doi:10.1038/nature13808')
    assert str(c) == s


# --- the RCSB Data API is opt-in and never fatal ---------------------------

RCSB_PAYLOAD = {
    'citation': [
        {'id': 'primary', 'rcsb_is_primary': 'Y',
         'title': 'Refined structure of human carbonic anhydrase II at 2.0 A resolution.',
         'rcsb_authors': ['Eriksson, A.E.', 'Jones, T.A.', 'Liljas, A.'],
         'rcsb_ORCID_identifiers': ['?', '?', '?'],
         'journal_abbrev': 'Proteins', 'journal_volume': '4',
         'page_first': '274', 'page_last': '282', 'year': 1988,
         'journal_id_ISSN': '0887-3585',
         'pdbx_database_id_DOI': '10.1002/prot.340040406',
         'pdbx_database_id_PubMed': 3151019},
        {'id': '1',
         'title': 'Crystallographic studies of inhibitor binding sites.',
         'rcsb_authors': ['Eriksson, A.E.', 'Kylsten, P.M.'],
         'journal_abbrev': 'Proteins', 'journal_volume': '4',
         'page_first': '283', 'year': 1988},
    ]
}


@pytest.fixture
def stub_rcsb(monkeypatch):
    """Answer the Data API from a canned payload; no socket is opened."""
    calls = []

    def fake_urlopen(url, timeout=None):
        calls.append(url)
        return io.BytesIO(json.dumps(RCSB_PAYLOAD).encode())

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    return calls


def test_citations_do_not_touch_the_network_by_default(monkeypatch):
    def explode(*a, **kw):
        raise AssertionError('citations() must not open a connection')

    monkeypatch.setattr(urllib.request, 'urlopen', explode)
    parse(PDB_1CA2).citations()
    parse(PDB_1CA2).citation_ids()
    parse(CIF_1CA2, 'mmCIF').citations()


def test_rcsb_citations_are_parsed(stub_rcsb):
    cits = citations_from_rcsb('1ca2')
    assert [c.role for c in cits] == [PRIMARY, RELATED]
    assert cits[0].last_page == '282'
    assert cits[0].orcids == []          # '?' placeholders dropped
    assert cits[0].source == 'rcsb-api'
    assert stub_rcsb == ['https://data.rcsb.org/rest/v1/core/entry/1ca2']


def test_enrichment_fills_gaps_without_discarding_local_data(stub_rcsb):
    p = PDBParser(source_db='rcsb', source_id='1ca2', filepath=PDB_1CA2).parse()
    plain = p.citations(role=PRIMARY)[0]
    rich = p.citations(role=PRIMARY, enrich=True)[0]
    assert plain.last_page == ''
    assert rich.last_page == '282'                              # gap filled
    assert rich.title.startswith('Refined structure')           # proper case
    assert rich.doi == plain.doi and rich.pmid == plain.pmid    # agreed values kept
    assert rich.source == 'JRNL+rcsb-api'


def test_a_failed_api_request_is_logged_not_raised(monkeypatch, caplog):
    def fail(*a, **kw):
        raise urllib.error.URLError('no route to host')

    monkeypatch.setattr(urllib.request, 'urlopen', fail)
    assert citations_from_rcsb('1ca2') == []
    p = PDBParser(source_db='rcsb', source_id='1ca2', filepath=PDB_1CA2).parse()
    assert p.citations(enrich=True) == p.citations()            # falls back intact


def test_merge_keeps_local_citations_a_remote_answer_omits():
    local = [Citation(id='primary', role=PRIMARY, doi='10.1/x'),
             Citation(id='7', role=RELATED, title='LOCAL ONLY')]
    remote = [Citation(id='primary', role=PRIMARY, last_page='9', source='rcsb-api')]
    merged = merge_citations(local, remote)
    assert len(merged) == 2
    assert merged[0].doi == '10.1/x' and merged[0].last_page == '9'
    assert merged[1].title == 'LOCAL ONLY'


def test_merge_with_nothing_remote_is_a_no_op():
    local = [Citation(id='primary', role=PRIMARY)]
    assert merge_citations(local, []) == local


# --- the two paths must agree on order and on 'no identifier' --------------

OUT_OF_ORDER_PAYLOAD = {
    'citation': [
        {'id': 'primary', 'rcsb_is_primary': 'Y', 'title': 'P'},
        {'id': '1', 'title': 'A'},
        {'id': '3', 'title': 'C'},
        {'id': '2', 'title': 'B'},
        {'id': '10', 'title': 'J'},
    ]
}


def test_numbered_citations_come_back_in_numeric_order(monkeypatch):
    """mmCIF does not oblige its citation rows to appear in id order -- 4HHB
    files them 1, 3, 4, 5, 6, 2, 7, 8 -- while REMARK 1 always counts up."""
    monkeypatch.setattr(urllib.request, 'urlopen',
                        lambda url, timeout=None: io.BytesIO(json.dumps(OUT_OF_ORDER_PAYLOAD).encode()))
    assert [c.id for c in citations_from_rcsb('xxxx')] == ['primary', '1', '2', '3', '10']


NO_PMID_PAYLOAD = {
    'citation': [{'id': 'primary', 'rcsb_is_primary': 'Y', 'title': 'P',
                  'pdbx_database_id_PubMed': -1}]
}


def test_the_minus_one_pubmed_sentinel_is_not_an_identifier(monkeypatch):
    """mmCIF and the RCSB API spell 'this paper has no PubMed ID' as -1 (1MBN,
    for one); a .pdb file simply carries no JRNL PMID line."""
    monkeypatch.setattr(urllib.request, 'urlopen',
                        lambda url, timeout=None: io.BytesIO(json.dumps(NO_PMID_PAYLOAD).encode()))
    assert citations_from_rcsb('xxxx')[0].pmid is None
