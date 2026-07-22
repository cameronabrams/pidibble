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


def test_unsupported_record_type_raises():
    """A type outside 1/3 is rejected rather than silently mis-emitted."""
    parser, w = _load()
    # COMPND is a type-2 (continuation) record
    try:
        w.emit(parser.parsed['COMPND'], 'COMPND')
    except NotImplementedError:
        return
    raise AssertionError('expected NotImplementedError for type-2 record')


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


def test_assembler_master_describes_output(tmp_path):
    """The regenerated MASTER counts match what was actually emitted."""
    parser, _ = _load()
    lines = parser.write_PDB(str(tmp_path / 'x.pdb'))
    master = next(l for l in lines if l.startswith('MASTER'))
    reparsed = PDBParser(filepath=str(tmp_path / 'x.pdb')).parse()
    m = reparsed.parsed['MASTER']
    # coordinate + geometry bookkeeping is reconstructable and must be right
    assert m.numCoord == len(parser.parsed['ATOM']) + len(parser.parsed['HETATM'])
    assert m.numTer == len(parser.parsed['TER'])
    assert m.numConect == len(parser.parsed['CONECT'])
    assert m.numXform == 6                       # ORIGX1-3 + SCALE1-3
    assert m.numHelix == len(parser.parsed['HELIX'])
    assert m.numSheet == len(parser.parsed['SHEET'])
    # REMARK/SEQRES are not yet writable, so the reduced file reports zero
    assert m.numRemark == 0 and m.numSeq == 0


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
