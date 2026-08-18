"""Round-trip tests for the fixed-column PDB writer (type-1/3 records).

The writer is the inverse of the parser, so the correctness gate is a
round-trip: parse a real entry, emit each record, and confirm the emitted line
(a) is byte-for-byte identical to the source line and (b) re-parses to the same
field values. Byte-exactness is the stronger of the two; the semantic check
guards records whose source spacing the parser already normalized away.
"""
import logging

import pytest

from pidibble.baseparsers import EmptyField
from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecord, PDBRecordList
from pidibble.pdbwrite import PDBWriter, FieldFormatter
from pidibble.hex import HexSerialEncoder, AtomSerialParser

logger = logging.getLogger(__name__)

# record types the prototype writer supports (all type-1 or type-3 in 4zmj)
ROUNDTRIP_KEYS = ['HEADER', 'CRYST1', 'ATOM', 'HETATM', 'SSBOND', 'LINK', 'CONECT']


def _load():
    # 4zmj.pdb is a committed fixture in this directory (pinned to the RCSB
    # release these byte-exact assertions were written against), so
    # PDBParser.fetch() finds it locally and never reaches the network.
    p = PDBParser(source_db='rcsb', source_id='4zmj').parse()
    w = PDBWriter(p.pdb_format_dict['record_formats'],
                  p.pdb_format_dict['custom_formats'])
    return p, w


def _records(parser, key):
    v = parser.parsed[key]
    return list(v) if isinstance(v, (list, PDBRecordList)) else [v]


def _assert_serials_equal(expected, got, label):
    """Compare two long serial-number sequences elementwise.

    Deliberately *not* a bare ``assert expected == got``: on the ~4500
    near-identical values of a full coordinate section, a mismatch sends
    pytest's assertion rewriting into difflib's ``_fancy_replace`` recursion,
    which takes effectively unbounded time and makes a plain regression look
    like a hung suite. Reporting a bounded slice fails just as loudly, in
    under a second, and points at the first bad serial.
    """
    assert len(got) == len(expected), f'{label}: got {len(got)} records, expected {len(expected)}'
    mismatch = [(i, e, g) for i, (e, g) in enumerate(zip(expected, got)) if e != g]
    assert not mismatch, f'{label}: {len(mismatch)} mismatches; first 10: {mismatch[:10]}'


def _source_lines():
    lines = {}
    with open('4zmj.pdb') as f:
        for line in f:
            lines.setdefault(line[:6].strip(), []).append(line.rstrip('\n').rstrip())
    return lines


def test_field_formatter_basics():
    f = FieldFormatter()
    assert f.format('', 'Float', 8, {}) == ' ' * 8
    assert f.format(-0.092, 'Float', 8, {}) == '  -0.092'          # default prec 3
    assert f.format(1.0, 'Float', 6, {'prec': 2}) == '  1.00'
    assert f.format(1, 'Integer', 5, {}) == '    1'                 # right by default
    assert f.format('LEU', 'String', 4, {}) == 'LEU '              # left by default
    assert f.format('N', 'String', 2, {'just': 'right'}) == ' N'
    # an absent typed field renders as blanks like a plain '' -- the writer must
    # not trip over the guard that makes it refuse numeric use
    assert f.format(EmptyField('length', 'Float'), 'Float', 5, {'prec': 2}) == ' ' * 5


def test_truncated_record_roundtrips_and_guards_its_empty_fields(tmp_path):
    """A record that stops short of its optional trailing fields (here SSBOND
    without sym1/sym2/length) parses, re-emits byte-identically, and hands back
    a guarded empty for the Float field rather than a '' that would multiply."""
    src = tmp_path / 'short.pdb'
    src.write_text('HEADER    TEST                                    01-JAN-00   XXXX\n'
                   'SSBOND   1 CYS A   54    CYS A   74\n'
                   'END\n')
    p = PDBParser(filepath=str(src)).parse()
    ssbond = p.parsed['SSBOND'][0]
    # the residue columns are canonical in this line, so they read correctly
    assert (ssbond.residue1.resName, ssbond.residue1.chainID, ssbond.residue1.seqNum) == ('CYS', 'A', 54)
    assert (ssbond.residue2.resName, ssbond.residue2.chainID, ssbond.residue2.seqNum) == ('CYS', 'A', 74)
    # a short record is not a nonconformance -- the columns are absent, not wrong
    assert p.nonconformances.types() == []
    # sym1/sym2 are String fields: plain '' is a legitimate value there
    assert ssbond.sym1 == '' and ssbond.sym2 == ''
    # length is declared Float, so it must refuse to behave like one
    assert ssbond.length == '' and not ssbond.length
    with pytest.raises(TypeError, match="'length'"):
        ssbond.length * 2
    lines = p.write_PDB(str(tmp_path / 'out.pdb'))
    assert [l for l in lines if l.startswith('SSBOND')] == ['SSBOND   1 CYS A   54    CYS A   74']


def test_atomname_justification():
    # 1-char element, short name -> indented one column (col 14)
    assert PDBWriter._emit_atomname('N', 'N', 4) == ' N  '
    assert PDBWriter._emit_atomname('CA', 'C', 4) == ' CA '
    # 2-char element -> starts in column 13
    assert PDBWriter._emit_atomname('FE', 'FE', 4) == 'FE  '
    # 4-char name fills the field
    assert PDBWriter._emit_atomname('HG21', 'H', 4) == 'HG21'


def test_roundtrip_byte_exact():
    """Every emitted line matches the corresponding source line exactly."""
    parser, w = _load()
    src = _source_lines()
    total = 0
    for key in ROUNDTRIP_KEYS:
        emitted = w.emit_all(_records(parser, key), key)
        srclines = src.get(key, [])
        assert len(emitted) == len(srclines), f'{key}: count {len(emitted)} != {len(srclines)}'
        for out, ref in zip(emitted, srclines):
            assert out == ref.rstrip(), f'{key}:\n  src={ref!r}\n  out={out!r}'
            total += 1
    logger.info(f'{total} records round-tripped byte-exact')
    assert total > 5000


def test_roundtrip_semantic():
    """Emitted lines re-parse to identical field values."""
    parser, w = _load()
    p = PDBParser(source_db='rcsb', source_id='4zmj')  # fresh mappers/format
    rf = p.pdb_format_dict['record_formats']
    checked = 0
    for key in ROUNDTRIP_KEYS:
        for rec in _records(parser, key):
            line = w.emit(rec, key)
            reparsed = PDBRecord.newrecord(key, line, rf[key], p.mappers)
            for field in rf[key]['fields']:
                a, b = getattr(rec, field), getattr(reparsed, field)
                if hasattr(a, '__dict__') and hasattr(b, '__dict__'):
                    assert a.__dict__ == b.__dict__, f'{key}.{field}: {a.__dict__} != {b.__dict__}'
                else:
                    assert a == b, f'{key}.{field}: {a!r} != {b!r}'
                checked += 1
    assert checked > 50000


def test_single_line_emit_rejects_multiline_record():
    """The single-line emit() refuses a continuation record (use emit_multiline);
    it must not silently mis-emit a multi-line record as one line."""
    parser, w = _load()
    try:
        w.emit(parser.parsed['COMPND'], 'COMPND')       # COMPND is type-2
    except NotImplementedError:
        return
    raise AssertionError('expected NotImplementedError from single-line emit()')


def test_multiline_emit_handles_continuation_record():
    """emit_multiline renders a type-2 record as one or more numbered lines."""
    parser, w = _load()
    out = w.emit_multiline(parser.parsed['KEYWDS'], 'KEYWDS')
    assert len(out) >= 1 and out[0].startswith('KEYWDS')
    assert 'HIV-1' in out[0]


# --- document assembler -----------------------------------------------------

# record types the assembler reconstructs and whose counts must survive a
# write / re-parse cycle
STRUCTURAL_KEYS = ['ATOM', 'HETATM', 'ANISOU', 'TER', 'CONECT',
                   'SSBOND', 'LINK', 'HELIX', 'SHEET', 'HET']


def test_assemble_roundtrip(tmp_path):
    """Parse -> write_PDB -> re-parse preserves the writable structure exactly."""
    parser, _ = _load()
    out = tmp_path / '4zmj_out.pdb'
    lines = parser.write_PDB(str(out))
    assert lines[-1] == 'END'
    assert any(l.startswith('MASTER') for l in lines)

    rewritten = PDBParser(filepath=str(out)).parse()
    for key in STRUCTURAL_KEYS:
        a = len(parser.parsed.get(key, []))
        b = len(rewritten.parsed.get(key, []))
        assert a == b, f'{key}: orig {a} != rewritten {b}'

    # coordinate values survive intact
    for i in (0, 2000, 4517):
        o, r = parser.parsed['ATOM'][i], rewritten.parsed['ATOM'][i]
        assert (o.x, o.y, o.z, o.occupancy, o.tempFactor) == \
               (r.x, r.y, r.z, r.occupancy, r.tempFactor)
        assert o.residue.__dict__ == r.residue.__dict__

    # crystallographic record survives intact
    ck = ['a', 'b', 'c', 'alpha', 'beta', 'gamma', 'sGroup', 'z']
    assert {k: parser.parsed['CRYST1'].__dict__[k] for k in ck} == \
           {k: rewritten.parsed['CRYST1'].__dict__[k] for k in ck}


def test_assembler_master_matches_original(tmp_path):
    """With types 2/4 written and REMARK/JRNL passed through, the regenerated
    MASTER equals the original entry's MASTER byte-for-byte."""
    parser, _ = _load()
    lines = parser.write_PDB(str(tmp_path / 'x.pdb'))
    master = next(l for l in lines if l.startswith('MASTER'))
    orig_master = next(l.rstrip() for l in open('4zmj.pdb') if l.startswith('MASTER'))
    assert master == orig_master

    reparsed = PDBParser(filepath=str(tmp_path / 'x.pdb')).parse()
    m = reparsed.parsed['MASTER']
    assert m.numCoord == len(parser.parsed['ATOM']) + len(parser.parsed['HETATM'])
    assert m.numTer == len(parser.parsed['TER'])
    assert m.numXform == 6                        # ORIGX1-3 + SCALE1-3
    assert m.numRemark == 649 and m.numSeq == 49  # passthrough + re-serialized


def test_continuation_and_group_records_reparse(tmp_path):
    """Type-2 (continuation) and type-4 (determinant-group) records survive a
    write / re-parse cycle with identical parsed content."""
    parser, _ = _load()
    parser.write_PDB(str(tmp_path / 'x.pdb'))
    q = PDBParser(filepath=str(tmp_path / 'x.pdb')).parse()

    # type-2: merged string / list / token content
    assert parser.parsed['TITLE'].title == q.parsed['TITLE'].title
    assert parser.parsed['KEYWDS'].keywds == q.parsed['KEYWDS'].keywds
    assert parser.parsed['COMPND'].get_token('CHAIN') == q.parsed['COMPND'].get_token('CHAIN')
    assert parser.parsed['SOURCE'].get_token('ORGANISM_TAXID') == \
           q.parsed['SOURCE'].get_token('ORGANISM_TAXID')

    # type-4: per-chain SEQRES residue lists
    for a, b in zip(parser.parsed['SEQRES'], q.parsed['SEQRES']):
        assert a.chainID == b.chainID and a.resNames == b.resNames

    # REVDAT: type-3 record that is genuinely multi-line (modNum 6 spans lines)
    orig6 = [r.records for r in parser.parsed['REVDAT'] if r.modNum == 6]
    new6 = [r.records for r in q.parsed['REVDAT'] if r.modNum == 6]
    assert orig6 == new6 and len(orig6[0]) == 6


def test_type4_seqres_byte_exact(tmp_path):
    """SEQRES re-serializes byte-exactly (13 residues/line matches the source)."""
    parser, _ = _load()
    lines = parser.write_PDB(str(tmp_path / 'x.pdb'))
    emitted = [l for l in lines if l.startswith('SEQRES')]
    src = [l.rstrip() for l in open('4zmj.pdb') if l.startswith('SEQRES')]
    assert emitted == src


def test_type6_passthrough_verbatim(tmp_path):
    """REMARK/JRNL are re-emitted verbatim from the source lines."""
    parser, _ = _load()
    lines = parser.write_PDB(str(tmp_path / 'x.pdb'))
    for card in ('REMARK', 'JRNL'):
        emitted = [l for l in lines if l[:6].strip() == card]
        src = [l.rstrip() for l in open('4zmj.pdb') if l[:6].strip() == card]
        assert emitted == src, f'{card} passthrough differs from source'


def test_full_key_preservation(tmp_path):
    """Every top-level parsed key survives a full write / re-parse round-trip."""
    parser, _ = _load()
    parser.write_PDB(str(tmp_path / 'x.pdb'))
    q = PDBParser(filepath=str(tmp_path / 'x.pdb')).parse()
    for key in parser.parsed:
        a, b = parser.parsed[key], q.parsed.get(key)
        assert b is not None, f'{key} missing after round-trip'
        la = len(a) if isinstance(a, (list, PDBRecordList)) else 1
        lb = len(b) if isinstance(b, (list, PDBRecordList)) else 1
        assert la == lb, f'{key}: orig {la} != rewritten {lb}'


def test_hex_serial_encoder_mirrors_parser():
    """The encoder is the exact inverse of AtomSerialParser, including the
    permanent decimal->hex trip and hex-encoded small back-references after it."""
    # emission order of a big file: serials climb past the ceiling, then a
    # CONECT section referencing a small atom (10) and a mid atom (255)
    serials = [99998, 99999, 100000, 100001, 250000, 1048575, 10, 255]
    enc, dec = HexSerialEncoder(), AtomSerialParser()
    for s in serials:
        assert dec(enc(s)) == s, f'{s} did not round-trip'
    # the boundary crossing is genuinely hex, not truncated decimal
    enc2 = HexSerialEncoder()
    assert enc2(99999) == '99999'
    assert enc2(100000) == '186A0'
    assert enc2(10) == 'A'            # tripped: small serial now hex


def test_big_serial_document_roundtrip(tmp_path):
    """A structure whose serials exceed 99999 writes hex serials that re-parse
    to the original decimal values across ATOM/HETATM/TER/CONECT."""
    parser, _ = _load()
    off = 100000
    for key in ('ATOM', 'HETATM', 'ANISOU', 'TER'):
        for r in parser.parsed[key]:
            r.serial += off
    for r in parser.parsed['CONECT']:
        r.serial += off
        for f in ('partner1', 'partner2', 'partner3', 'partner4'):
            v = getattr(r, f)
            if v != '':
                setattr(r, f, v + off)

    out = tmp_path / 'big.pdb'
    lines = parser.write_PDB(str(out))
    # the first atom's serial is emitted as hex, not clipped decimal
    assert next(l for l in lines if l.startswith('ATOM'))[6:11] == '186A1'

    q = PDBParser(filepath=str(out)).parse()
    for key in ('ATOM', 'TER'):
        _assert_serials_equal([r.serial for r in parser.parsed[key]],
                              [r.serial for r in q.parsed[key]], key)

    def conect(rec):
        return (rec.serial,) + tuple(getattr(rec, f)
                                     for f in ('partner1', 'partner2', 'partner3', 'partner4'))
    _assert_serials_equal(sorted(conect(r) for r in parser.parsed['CONECT']),
                          sorted(conect(r) for r in q.parsed['CONECT']), 'CONECT')


def test_assembler_ter_cards(tmp_path):
    """TER cards are reconstructed at chain breaks with correct residues."""
    parser, _ = _load()
    lines = parser.write_PDB(str(tmp_path / 'x.pdb'))
    ters = [l for l in lines if l.startswith('TER')]
    assert ters == ['TER    3544      VAL G 505', 'TER    4520      ASP B 664']


def test_assembler_without_anisou(tmp_path):
    """anisou=False drops the anisotropic records from the coordinate section."""
    parser, _ = _load()
    lines = parser.write_PDB(str(tmp_path / 'x.pdb'), anisou=False)
    assert not any(l.startswith('ANISOU') for l in lines)
    assert sum(l.startswith('ATOM') for l in lines) == len(parser.parsed['ATOM'])


# --- CHARMM write/read dialect ----------------------------------------------

CHARMM_FIXTURE = 'charmm_glycan.pdb'   # wide resNames + segID + iCode + altLoc + TER

# fields pestifer relies on surviving the CHARMM round-trip
CHARMM_FIELDS = ('serial', 'name', 'altLoc', 'x', 'y', 'z',
                 'occupancy', 'tempFactor', 'segID', 'element', 'charge')


def test_charmm_dialect_parses_wide_resname_and_segid():
    """The CHARMM dialect reads 6-char resNames, the segID column, and
    leading-digit atom names that the standard 3-4 column model can't."""
    p = PDBParser(filepath=CHARMM_FIXTURE, dialect='charmm').parse()
    a = p.parsed['ATOM'][0]
    assert a.residue.resName == 'BGLCNA'          # 6 chars, not truncated
    assert a.residue.chainID == 'A' and a.residue.seqNum == 1
    assert a.segID == 'GLYC'
    assert (a.x, a.y, a.z) == (12.345, 23.456, 34.567)
    assert p.parsed['ATOM'][1].name == '1HD1'     # leading-digit 4-char name
    assert p.parsed['ATOM'][1].residue.resName == 'ANE5AC'
    assert p.parsed['ATOM'][2].altLoc == 'A'      # altLoc preserved
    assert p.parsed['ATOM'][4].residue.iCode == 'A'  # insertion code preserved


def test_charmm_dialect_roundtrip(tmp_path):
    """Acceptance test: parse -> write -> parse is the identity for every field
    pestifer relies on, on a CHARMM structure with glycans and a segID column."""
    p = PDBParser(filepath=CHARMM_FIXTURE, dialect='charmm').parse()
    out = tmp_path / 'charmm_out.pdb'
    p.write_PDB(str(out), anisou=False, dialect='charmm')
    q = PDBParser(filepath=str(out), dialect='charmm').parse()

    for key in ('ATOM', 'HETATM'):
        assert len(p.parsed[key]) == len(q.parsed[key])
        for a, b in zip(p.parsed[key], q.parsed[key]):
            for f in CHARMM_FIELDS:
                assert getattr(a, f) == getattr(b, f), \
                    f'{key}.{f}: {getattr(a, f)!r} != {getattr(b, f)!r}'
            assert a.residue.__dict__ == b.residue.__dict__   # resName/chain/seq/iCode
    # TER cards (with wide resNames) survive too
    assert len(p.parsed['TER']) == len(q.parsed['TER']) == 2
    assert [t.residue.resName for t in q.parsed['TER']] == ['AGLC', 'HSD']


def test_charmm_coordinate_columns_pinned(tmp_path):
    """The hazard is designed out: x/y/z stay at columns 31-54 and segID at
    73-76 regardless of resName width, so fixed-column readers never scramble."""
    p = PDBParser(filepath=CHARMM_FIXTURE, dialect='charmm').parse()
    lines = p.write_PDB(str(tmp_path / 'out.pdb'), anisou=False, dialect='charmm')
    coord = [l for l in lines if l[:6].strip() in ('ATOM', 'HETATM')]
    for l in coord:
        padded = l.ljust(80)
        assert float(padded[30:38]) == float(padded[30:38])   # x parses at 31-38
        # x column is numeric and the resName never bleeds into it
        assert padded[30:54].strip() and padded[24:30].strip() != ''
    first = next(l for l in coord if l.startswith('ATOM')).ljust(80)
    assert first[17:23] == 'BGLCNA'      # resName at 18-23
    assert first[30:38] == '  12.345'    # x pinned at 31-38
    assert first[72:76] == 'GLYC'        # segID at 73-76


def test_dialect_selects_format_table():
    """Constructing with a dialect swaps the coordinate-record column model;
    'standard' has no segID field, 'charmm' does and uses the wide residue."""
    std = PDBParser(dialect='standard')
    chm = PDBParser(dialect='charmm')
    assert 'segID' not in std.record_formats['ATOM']['fields']
    assert std.record_formats['ATOM']['fields']['residue'][0] == 'Residue10'
    assert 'segID' in chm.record_formats['ATOM']['fields']
    assert chm.record_formats['ATOM']['fields']['residue'][0] == 'ResidueCharmm'
    # x stays pinned at column 31 in both dialects
    assert std.record_formats['ATOM']['fields']['x'][1] == [31, 38]
    assert chm.record_formats['ATOM']['fields']['x'][1] == [31, 38]


def test_unknown_dialect_rejected():
    try:
        PDBParser(dialect='amber')
    except ValueError:
        return
    raise AssertionError('expected ValueError for an unknown dialect')


# --- JRNL serialized from mmCIF -------------------------------------------
#
# A .pdb source round-trips JRNL by passthrough (test_type6_passthrough_verbatim
# above). An mmCIF source has no JRNL lines to pass through, so the block is
# serialized from the JRNL.* sub-records the `citation` category supplies. The
# correctness gate is the entry's own .pdb file: the two are committed
# fixtures of the same release, so the emitted block should match it exactly.
# They live in ../test_rcsb rather than being duplicated here (1ca2.cif alone
# is 276K).

CIF_1CA2 = '../test_rcsb/1ca2.cif'
PDB_1CA2 = '../test_rcsb/1ca2.pdb'
CIF_4ZMJ = '../test_rcsb/4zmj.cif'
PDB_4ZMJ = '../test_rcsb/4zmj.pdb'


def _jrnl(lines):
    return [l for l in lines if l[:6].strip() == 'JRNL']


def _source_jrnl(path):
    return [l.rstrip() for l in open(path) if l[:6].strip() == 'JRNL']


def _emit_jrnl_block(parser):
    """Emit just the JRNL block, without assembling a whole document.

    Full CIF-to-PDB assembly is unrelatedly broken for some entries (mmCIF puts
    struct_conn.id, e.g. 'disulf1', into the Integer SSBOND.serNum), so the
    tests that need only the JRNL block drive the writer directly.
    """
    from pidibble.pdbwrite import JRNL_SUBKEYS
    w = PDBWriter(parser.record_formats, parser.pdb_format_dict['custom_formats'])
    subfmts = parser.record_formats['JRNL']['subrecords']['formats']
    out = []
    for subkey in JRNL_SUBKEYS:
        rec = parser.parsed.get(f'JRNL.{subkey}')
        if rec is not None and subkey in subfmts:
            out.extend(w.emit_subrecord(rec, 'JRNL', subkey, subfmts[subkey]))
    return out


def test_jrnl_from_mmcif_is_byte_identical_to_the_pdb_file(tmp_path):
    """The whole JRNL block, rebuilt from mmCIF, matches the entry's own .pdb."""
    p = PDBParser(input_format='mmCIF', filepath=CIF_1CA2).parse()
    emitted = _jrnl(p.write_PDB(str(tmp_path / 'x.pdb')))
    assert emitted == _source_jrnl(PDB_1CA2)
    assert len(emitted) == 7          # AUTH, TITL x2, REF, REFN, PMID, DOI


def test_jrnl_from_mmcif_wraps_a_long_author_list_like_the_pdb():
    """4ZMJ has 53 authors across 10 AUTH lines. The separator trails the line
    it ends, and packing keeps a column free for it."""
    p = PDBParser(input_format='mmCIF', filepath=CIF_4ZMJ).parse()
    emitted = _emit_jrnl_block(p)
    want = _source_jrnl(PDB_4ZMJ)

    auth = [l for l in emitted if l[12:16] == 'AUTH']
    assert len(auth) == 10
    assert [l for l in want if l[12:16] == 'AUTH'] == auth
    # every continued line ends with the separator; the last one does not
    assert all(l.rstrip().endswith(',') for l in auth[:-1])
    assert not auth[-1].rstrip().endswith(',')
    assert all(len(l.rstrip()) <= 80 for l in emitted)

    # the entry's ESSN is the sole thing mmCIF cannot supply (it records an
    # ISSN without saying whether it is the print or the electronic one)
    diffs = [(a, b) for a, b in zip(emitted, want) if a != b]
    assert len(emitted) == len(want)
    assert diffs == [(next(l for l in emitted if l[12:16] == 'REFN'),
                      next(l for l in want if l[12:16] == 'REFN'))]
    assert diffs[0][0].replace('ISSN', 'ESSN') == diffs[0][1]


def test_jrnl_from_mmcif_reparses_to_the_same_records(tmp_path):
    """The emitted block re-parses to the sub-records it was built from."""
    p = PDBParser(input_format='mmCIF', filepath=CIF_1CA2).parse()
    # re-parse the emitted JRNL block on its own: a whole assembled document
    # would trip over unrelated CIF-to-PDB breakage (see _emit_jrnl_block)
    out = tmp_path / 'jrnl.pdb'
    out.write_text('\n'.join(_emit_jrnl_block(p)) + '\nEND\n')
    q = PDBParser(filepath=str(out)).parse()
    assert q.parsed['JRNL.AUTH'].authorList == p.parsed['JRNL.AUTH'].authorList
    assert q.parsed['JRNL.TITL'].title == p.parsed['JRNL.TITL'].title
    assert q.parsed['JRNL.REF'].year == p.parsed['JRNL.REF'].year
    assert q.parsed['JRNL.REF'].volume == p.parsed['JRNL.REF'].volume
    assert q.parsed['JRNL.PMID'].pmid == p.parsed['JRNL.PMID'].pmid
    assert q.parsed['JRNL.DOI'].doi == p.parsed['JRNL.DOI'].doi
    # the citation survives the write cycle in everything the PDB format can
    # carry. Title case and last_page are not among those things -- the mmCIF
    # holds them and a .pdb cannot, which is the documented asymmetry.
    before, after = p.citations()[0], q.citations()[0]
    assert (after.authors, after.doi, after.pmid, after.year) == \
           (before.authors, before.doi, before.pmid, before.year)
    assert (after.journal, after.volume, after.first_page) == \
           (before.journal.upper(), before.volume, before.first_page)
    assert after.title == before.title.upper()


def test_pdb_source_still_passes_jrnl_through_verbatim(tmp_path):
    """Serialization is the mmCIF fallback only: a .pdb source keeps the
    byte-exact passthrough, even where the two would differ."""
    p = PDBParser(filepath=PDB_4ZMJ).parse()
    emitted = _jrnl(p.write_PDB(str(tmp_path / 'x.pdb')))
    assert emitted == _source_jrnl(PDB_4ZMJ)
    assert any('ESSN' in l for l in emitted)      # the passthrough kept it


def test_remark_is_still_omitted_for_an_mmcif_source(tmp_path, caplog):
    """REMARK is not regenerated: its free text is held as an LList whose line
    breaks the parser discards, so it cannot be reproduced faithfully."""
    p = PDBParser(input_format='mmCIF', filepath=CIF_1CA2).parse()
    with caplog.at_level(logging.INFO):
        lines = p.write_PDB(str(tmp_path / 'x.pdb'))
    assert not [l for l in lines if l[:6].strip() == 'REMARK']
    assert _jrnl(lines)                            # but JRNL is there
    assert any('REMARK' in r.message for r in caplog.records)


def test_break_long_token_prefers_a_comma_boundary():
    """A journal name too wide for its column breaks after a comma rather than
    being clipped: PHILOS.TRANS.R.SOC.LONDON, / SER.B (4INS)."""
    assert PDBWriter._break_long_token('PHILOS.TRANS.R.SOC.LONDON,SER.B', 28) == \
           ['PHILOS.TRANS.R.SOC.LONDON,', 'SER.B']
    # no comma to break at: a hard break, but nothing is lost
    assert ''.join(PDBWriter._break_long_token('A' * 70, 28)) == 'A' * 70
    assert PDBWriter._break_long_token('SHORT', 28) == ['SHORT']


def test_pack_items_reserves_a_column_only_while_items_remain():
    """The trailing separator is owed on every line but the last, so the last
    line may use the full width."""
    # 'aaaa,bbbb' is 9 wide; with a trailing ',' it would be 10
    assert PDBWriter._pack_items(['aaaa', 'bbbb'], ',', 9, reserve=1) == [['aaaa', 'bbbb']]
    # a third item makes the first line owe its comma, forcing the wrap earlier
    assert PDBWriter._pack_items(['aaaa', 'bbbb', 'cc'], ',', 9, reserve=1) == \
           [['aaaa'], ['bbbb', 'cc']]


# --- a whole document written from mmCIF ----------------------------------

def test_cif_sourced_document_round_trips(tmp_path):
    """mmCIF -> write_PDB -> re-parse. Previously this raised: the mmCIF
    mapping put struct_conn.id ('disulf1') into the Integer SSBOND.serNum."""
    p = PDBParser(input_format='mmCIF', filepath=CIF_4ZMJ).parse()
    q = PDBParser(filepath=str(tmp_path / 'x.pdb')).parse() if \
        p.write_PDB(str(tmp_path / 'x.pdb')) else None
    assert q is not None
    assert [r.serNum for r in p.parsed['SSBOND']] == [r.serNum for r in q.parsed['SSBOND']]
    # the entry's paper survives the whole trip
    assert [(c.doi, c.pmid) for c in q.citation_ids(role='primary')] == \
           [(c.doi, c.pmid) for c in p.citation_ids(role='primary')]


def test_contentless_card_is_omitted_not_written_blank(tmp_path, caplog):
    """An mmCIF COMPND carries molID/molName/chains, not the raw `compound`
    SList the PDB declares, so there is nothing to write. Emitting a bare
    'COMPND' produced a non-conformant record that broke re-parsing."""
    p = PDBParser(input_format='mmCIF', filepath=CIF_1CA2).parse()
    assert 'COMPND' in p.parsed                      # the record is there
    with caplog.at_level(logging.INFO):
        lines = p.write_PDB(str(tmp_path / 'x.pdb'))
    for card in ('COMPND', 'SOURCE'):
        assert not [l for l in lines if l[:6].strip() == card], f'blank {card} written'
    assert any('no writable content' in r.message for r in caplog.records)
    PDBParser(filepath=str(tmp_path / 'x.pdb')).parse()   # re-parses cleanly


def test_a_record_with_no_continued_content_but_other_fields_is_kept(tmp_path):
    """The omission rule tests the rendered line, not one field: a REVDAT
    initial-release entry has an empty `records` list but real modNum/date."""
    parser, _ = _load()
    lines = parser.write_PDB(str(tmp_path / 'x.pdb'))
    q = PDBParser(filepath=str(tmp_path / 'x.pdb')).parse()
    assert len(q.parsed['REVDAT']) == len(parser.parsed['REVDAT'])
    initial = [r for r in q.parsed['REVDAT'] if not r.records]
    assert initial and all(r.modNum and r.modDate for r in initial)


def test_parse_tokens_survives_a_non_list_field(caplog):
    """A malformed tokenized card must not abort the parse. It used to assert."""
    fmt = {'fields': {'continuation': ['String', [9, 10]],
                      'compound': ['SList', [11, 80]]},
           'token_formats': {'compound': {'tokens': {'MOL_ID': {'type': 'Integer'}}}}}
    rec = PDBRecord({'key': 'COMPND', 'compound': '', 'format': fmt})
    with caplog.at_level(logging.WARNING):
        rec.parse_tokens({})
    assert rec.tokengroups == {'compound': {}}
    assert any('not a list of token-strings' in r.message for r in caplog.records)
