"""Round-trip tests for the fixed-column PDB writer (type-1/3 records).

The writer is the inverse of the parser, so the correctness gate is a
round-trip: parse a real entry, emit each record, and confirm the emitted line
(a) is byte-for-byte identical to the source line and (b) re-parses to the same
field values. Byte-exactness is the stronger of the two; the semantic check
guards records whose source spacing the parser already normalized away.
"""
import logging

from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecord, PDBRecordList
from pidibble.pdbwrite import PDBWriter, FieldFormatter
from pidibble.hex import HexSerialEncoder, AtomSerialParser

logger = logging.getLogger(__name__)

# record types the prototype writer supports (all type-1 or type-3 in 4zmj)
ROUNDTRIP_KEYS = ['HEADER', 'CRYST1', 'ATOM', 'HETATM', 'SSBOND', 'LINK', 'CONECT']


def _load():
    p = PDBParser(source_db='rcsb', source_id='4zmj').parse()
    w = PDBWriter(p.pdb_format_dict['record_formats'],
                  p.pdb_format_dict['custom_formats'])
    return p, w


def _records(parser, key):
    v = parser.parsed[key]
    return list(v) if isinstance(v, (list, PDBRecordList)) else [v]


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
    assert [r.serial for r in parser.parsed['ATOM']] == [r.serial for r in q.parsed['ATOM']]
    assert [r.serial for r in parser.parsed['TER']] == [r.serial for r in q.parsed['TER']]

    def conect(rec):
        return (rec.serial,) + tuple(getattr(rec, f)
                                     for f in ('partner1', 'partner2', 'partner3', 'partner4'))
    assert sorted(conect(r) for r in parser.parsed['CONECT']) == \
           sorted(conect(r) for r in q.parsed['CONECT'])


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
