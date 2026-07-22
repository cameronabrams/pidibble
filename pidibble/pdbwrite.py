# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""

.. module:: pdbwrite
   :synopsis: fixed-column writer for parsed PDB records (type-1/3 prototype)

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

This is the inverse of :class:`.baseparsers.StringParser`: it turns parsed
record objects back into spec-conformant, fixed-column PDB lines.

Parsing is driven by ``[typestring, [start, end]]`` field specs; writing is
driven by the *same* specs plus an optional third element carrying formatting
hints the parser discards but the writer needs::

    x:          [Float, [31, 38]]                 # default: 3 decimals, right
    occupancy:  [Float, [55, 60], {prec: 2}]      # 2 decimals
    name:       [String, [13, 16], {just: atomname}]
    element:    [String, [77, 78], {just: right}]

Scope of this prototype: type-1 (single-line) and type-3 (multiple single-line)
records, including composite fields (``Residue10`` &c.). Continuation records
(type 2), grouped/embedded records (types 4-6), and hex serials > 99999 are out
of scope and raise or degrade gracefully.
"""
import logging

from .pdbrecord import PDBRecordList

logger = logging.getLogger(__name__)

#: default number of decimal places for a Float field with no ``prec`` hint
DEFAULT_FLOAT_PREC = 3

#: type-5 (grouping) records that are nonetheless a single emittable line
TYPE5_EMITTABLE = {'TER'}


def _hints(spec):
    """Return the optional formatting-hints dict from a field spec (or {})."""
    return spec[2] if len(spec) > 2 else {}


class FieldFormatter:
    """
    Inverse of the parse-time type coercion: value -> fixed-width field string.

    A formatter is the writer-side counterpart of the ``typemap`` used by
    :class:`.baseparsers.StringParser`. Given a value, its declared type, the
    field width (from the byte range), and any hints, it returns a string of
    *exactly* ``width`` characters.
    """

    def format(self, value, typestring, width, hints):
        """
        Format one scalar field value into exactly ``width`` characters.

        Parameters
        ----------
        value : any
            The parsed value (``''`` or ``None`` means an empty field).
        typestring : str
            The field's declared type (``'String'``, ``'Integer'``,
            ``'HxInteger'``, ``'Float'``).
        width : int
            The field width in columns.
        hints : dict
            Optional ``{'prec': int, 'just': 'left'|'right'}`` formatting hints.

        Returns
        -------
        str
            A string of length exactly ``width``.
        """
        just = hints.get('just')
        if value == '' or value is None:
            return ' ' * width

        if typestring == 'Float':
            prec = hints.get('prec', DEFAULT_FLOAT_PREC)
            s = f'{float(value):.{prec}f}'
            return self._fit(s, width, just or 'right')

        if typestring in ('Integer', 'HxInteger'):
            iv = int(value)
            if typestring == 'HxInteger' and iv > 99999:
                # hex-serial writing (files > 99999 atoms) is out of scope here
                logger.warning(f'serial {iv} exceeds 99999; hex writing not implemented')
            return self._fit(str(iv), width, just or 'right')

        # String
        return self._fit(str(value), width, just or 'left')

    @staticmethod
    def _fit(s, width, just):
        """Pad or (on overflow) clip ``s`` to exactly ``width`` for the given justification."""
        if len(s) > width:
            # keep the significant end: rightmost for right-justified fields
            return s[-width:] if just == 'right' else s[:width]
        return s.rjust(width) if just == 'right' else s.ljust(width)


class PDBWriter:
    """
    Emit spec-conformant fixed-column PDB lines from parsed records.

    Parameters
    ----------
    format_dict : dict
        The ``record_formats`` mapping (as loaded from ``pdb_format.yaml``).
    custom_formats : dict
        The ``custom_formats`` mapping (composite field layouts).
    """

    def __init__(self, format_dict, custom_formats):
        self.record_formats = format_dict
        self.custom_formats = custom_formats
        self.formatter = FieldFormatter()

    def emit(self, record, key=None):
        """
        Render one parsed record as a single 80-column line (right-trimmed).

        Parameters
        ----------
        record : BaseRecord
            A parsed record whose attributes are the fields named in its format.
        key : str, optional
            The record key (e.g. ``'ATOM'``). Defaults to ``record.key``; a
            dotted subrecord key (``'REMARK.350'``) uses its head as the card
            name.

        Returns
        -------
        str
            The formatted record line, trailing spaces removed.
        """
        key = key or getattr(record, 'key', None)
        card = key.split('.')[0]
        fmt = self.record_formats[card]
        rtype = fmt.get('type')
        if rtype not in (1, 3) and card not in TYPE5_EMITTABLE:
            raise NotImplementedError(
                f'{card} is record type {rtype}; only type-1/3 writing is implemented')

        line = bytearray(b' ' * 80)
        self._place(line, 1, 6, card.ljust(6))

        for fname, spec in fmt['fields'].items():
            typestring, (start, end) = spec[0], spec[1]
            hints = _hints(spec)
            width = end - start + 1
            value = getattr(record, fname, '')

            if typestring in self.custom_formats:
                text = self._emit_composite(typestring, value, width)
            elif hints.get('just') == 'atomname':
                text = self._emit_atomname(value, getattr(record, 'element', ''), width)
            else:
                text = self.formatter.format(value, typestring, width, hints)

            self._place(line, start, end, text)

        return line.decode('ascii').rstrip()

    def emit_all(self, records, key):
        """Render a :class:`.pdbrecord.PDBRecordList` (or list) of records."""
        return [self.emit(r, key) for r in records]

    # ---- internals -------------------------------------------------------

    def _emit_composite(self, typestring, subrecord, width):
        """Assemble a composite field (e.g. ``Residue10``) from its sub-record."""
        subfmt = self.custom_formats[typestring]
        buf = bytearray(b' ' * width)
        for sfname, sspec in subfmt.items():
            styp, (s, e) = sspec[0], sspec[1]
            sval = getattr(subrecord, sfname, '')
            text = self.formatter.format(sval, styp, e - s + 1, _hints(sspec))
            self._place(buf, s, e, text)
        return buf.decode('ascii')

    @staticmethod
    def _emit_atomname(name, element, width):
        """
        Apply the PDB atom-name justification rule within cols 13-16.

        Two-character element symbols (FE, CL, ...) and four-character atom
        names start in column 13; otherwise the name is indented one column so
        the (one-character) element aligns under column 14.
        """
        name = str(name)
        if len(name) >= width or len(str(element)) == 2:
            return name.ljust(width)[:width]
        return (' ' + name).ljust(width)[:width]

    @staticmethod
    def _place(buf, start, end, text):
        """Drop ``text`` (already sized) into 1-based columns ``[start, end]``."""
        width = end - start + 1
        assert len(text) == width, f'field text {text!r} is {len(text)} chars, expected {width}'
        buf[start - 1:end] = text.encode('ascii')


#: emission order for the sections preceding the coordinate block. Only the
#: type-1/3 members are actually written; the rest are listed so the order is
#: legible and future-proof, and are silently skipped when not emittable.
_PRE_COORD_ORDER = [
    'HEADER', 'OBSLTE', 'TITLE', 'SPLIT', 'CAVEAT', 'COMPND', 'SOURCE', 'KEYWDS',
    'EXPDTA', 'NUMMDL', 'MDLTYP', 'AUTHOR', 'REVDAT', 'SPRSDE', 'JRNL',
    'DBREF', 'DBREF1', 'DBREF2', 'SEQADV', 'SEQRES', 'MODRES',
    'HET', 'HETNAM', 'HETSYN', 'FORMUL',
    'HELIX', 'SHEET',
    'SSBOND', 'LINK', 'CISPEP',
    'SITE',
    'CRYST1', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3',
    'MTRIX1', 'MTRIX2', 'MTRIX3',
]

#: records whose emitted counts feed the MASTER bookkeeping line
_XFORM_KEYS = ['ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3',
               'MTRIX1', 'MTRIX2', 'MTRIX3']


def assemble_pdb(parser, anisou=True, include_master=True):
    """
    Assemble a conformant PDB document from a parsed structure.

    Emits every type-1/3 record (plus ``TER``) present in ``parser.parsed`` in
    canonical section order, reconstructs the coordinate section (``ATOM`` with
    interleaved ``ANISOU`` and chain-terminating ``TER`` cards, then
    ``HETATM``), and regenerates the ``MASTER`` and ``END`` bookkeeping records
    from the emitted content. Record types that cannot yet be written
    (continuation, grouped, embedded — e.g. ``REMARK``, ``COMPND``, ``SEQRES``)
    are skipped and reported via the logger, so the result is a *reduced* but
    internally consistent file.

    Parameters
    ----------
    parser : PDBParser
        A parsed structure (``parser.parsed`` populated).
    anisou : bool, optional
        Interleave ``ANISOU`` records after their matching atom (default True).
    include_master : bool, optional
        Regenerate a ``MASTER`` record from the emitted counts (default True).

    Returns
    -------
    list of str
        The document as a list of record lines (no trailing newlines).
    """
    rf = parser.pdb_format_dict['record_formats']
    cf = parser.pdb_format_dict['custom_formats']
    w = PDBWriter(rf, cf)
    parsed = parser.parsed

    if 'MODEL' in parsed:
        logger.warning('multi-model coordinate sections are not reconstructed; '
                       'emitting atoms as a single flat model')

    lines = []
    counts = {}

    def _emit_key(key):
        recs = parsed.get(key)
        if recs is None:
            return
        if rf.get(key, {}).get('type') not in (1, 3):
            return
        rlist = list(recs) if isinstance(recs, (list, PDBRecordList)) else [recs]
        for r in rlist:
            lines.append(w.emit(r, key))
        counts[key] = counts.get(key, 0) + len(rlist)

    # --- everything up to and including the crystallographic section --------
    for key in _PRE_COORD_ORDER:
        _emit_key(key)

    # --- coordinate section -------------------------------------------------
    anisou_map = {}
    if anisou:
        for a in parsed.get('ANISOU', []):
            anisou_map[a.serial] = a
    ter_by_chain = {t.residue.chainID: t for t in parsed.get('TER', [])}

    def _emit_atom(rec, key):
        lines.append(w.emit(rec, key))
        counts[key] = counts.get(key, 0) + 1
        anis = anisou_map.get(rec.serial)
        if anis is not None:
            lines.append(w.emit(anis, 'ANISOU'))

    atoms = list(parsed.get('ATOM', []))
    for i, a in enumerate(atoms):
        _emit_atom(a, 'ATOM')
        next_chain = atoms[i + 1].residue.chainID if i + 1 < len(atoms) else None
        if a.residue.chainID != next_chain:               # chain boundary
            ter = ter_by_chain.get(a.residue.chainID)
            if ter is not None:
                lines.append(w.emit(ter, 'TER'))
                counts['TER'] = counts.get('TER', 0) + 1

    for h in parsed.get('HETATM', []):
        _emit_atom(h, 'HETATM')

    # --- connectivity and bookkeeping --------------------------------------
    _emit_key('CONECT')
    if include_master:
        lines.append(_master_line(w, parsed, counts))
    lines.append('END')

    _report_skipped(parser, counts)
    return lines


def _master_line(writer, parsed, counts):
    """Build the MASTER record describing the emitted content."""
    numxform = sum(counts.get(k, 0) for k in _XFORM_KEYS)
    fields = {
        'numRemark': counts.get('REMARK', 0),
        'numHet': counts.get('HET', 0),
        'numHelix': counts.get('HELIX', 0),
        'numSheet': counts.get('SHEET', 0),
        'numTurn': 0,
        'numSite': counts.get('SITE', 0),
        'numXform': numxform,
        'numCoord': counts.get('ATOM', 0) + counts.get('HETATM', 0),
        'numTer': counts.get('TER', 0),
        'numConect': counts.get('CONECT', 0),
        'numSeq': counts.get('SEQRES', 0),
    }
    rec = _SimpleRecord(fields)
    line = writer.emit(rec, 'MASTER')
    # cols 16-20 are a spec-mandated "0" that pidibble's field map omits
    line = f'{line[:15]}{"0":>5}{line[20:]}'.rstrip()
    return line


def _report_skipped(parser, counts):
    """Log the record types present in the parse but absent from the output."""
    present = {}
    for key in parser.parsed:
        card = key.split('.')[0]
        present[card] = present.get(card, 0) + 1
    skipped = sorted(c for c in present if c not in counts
                     and c not in ('MASTER', 'END', 'MODEL', 'ENDMDL', 'ANISOU'))
    if skipped:
        logger.info(f'assemble_pdb skipped {len(skipped)} non-writable record '
                    f'type(s) (type 2/4/5/6): {", ".join(skipped)}')


class _SimpleRecord:
    """A minimal attribute bag standing in for a parsed record (for MASTER)."""

    def __init__(self, fields):
        self.__dict__.update(fields)
