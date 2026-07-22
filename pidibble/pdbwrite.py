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

from .hex import HexSerialEncoder
from .pdbrecord import PDBRecordList

logger = logging.getLogger(__name__)

#: default number of decimal places for a Float field with no ``prec`` hint
DEFAULT_FLOAT_PREC = 3

#: type-5 (grouping) records that are nonetheless a single emittable line
TYPE5_EMITTABLE = {'TER'}

#: how a list-valued field is re-joined into its column, by list type
LIST_SEP = {'CList': ', ', 'SList': '; ', 'WList': ' ', 'DList': ': ', 'LList': ' '}


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
            # stateless fallback: hex once past the decimal ceiling (the
            # document path uses a stateful encoder instead — see PDBWriter)
            s = format(iv, 'X') if (typestring == 'HxInteger' and iv > 99999) else str(iv)
            return self._fit(s, width, just or 'right')

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
        # stateful decimal->hex serial encoder, threaded through every serial
        # field this writer emits in order (see HexSerialEncoder)
        self.hex = HexSerialEncoder()

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

        values = {f: getattr(record, f, '') for f in fmt['fields']}
        return self._render_line(card, values, fmt)

    def emit_all(self, records, key):
        """Render a :class:`.pdbrecord.PDBRecordList` (or list) of records."""
        return [self.emit(r, key) for r in records]

    def emit_multiline(self, record, key=None):
        """
        Render a continuation (type-2) or determinant-group (type-4) record as
        one or more numbered lines.

        Both record families share a shape: a *counter* field, some *constant*
        fields repeated on every line, and a single *continued* field whose
        content is chunked across the lines. The differences are all derivable
        from the format: whether the counter is the blank-first ``continuation``
        field or a one-based serial (``serNum``/``seqNum``), and whether the
        continued field is a real column or a ``concatenate`` of sub-fields
        (e.g. ``SITE``'s four residue slots).

        Parameters
        ----------
        record : PDBRecord
            A parsed (already-merged) record.
        key : str, optional
            The record key; defaults to ``record.key``.

        Returns
        -------
        list of str
            The record rendered as one line per chunk.
        """
        key = key or getattr(record, 'key', None)
        card = key.split('.')[0]
        fmt = self.record_formats[card]
        continues = fmt.get('continues', [])
        assert len(continues) == 1, f'{card}: expected exactly one continued field'
        cf = continues[0]
        counter = fmt.get('continuation', 'continuation')
        blank_first = counter == 'continuation'
        concat = fmt.get('concatenate', {})

        # a String counter (the blank-first `continuation` field) must be
        # right-justified by hand so its number aligns like the spec expects;
        # an Integer counter (serNum/seqNum) is right-justified by the formatter
        counter_type, (cs, ce) = fmt['fields'][counter][0], fmt['fields'][counter][1]
        counter_w = ce - cs + 1

        chunks = self._chunk(cf, getattr(record, cf, ''), fmt, concat)

        # constant fields = everything that is neither the counter nor part of
        # the continued content
        content = set(concat.get(cf, [cf])) | {cf}
        const = {f: getattr(record, f, '') for f in fmt['fields']
                 if f != counter and f not in content}

        lines = []
        for i, chunk in enumerate(chunks):
            values = dict(const)
            if blank_first and i == 0:
                values[counter] = ''
            else:
                values[counter] = str(i + 1).rjust(counter_w) \
                    if counter_type == 'String' else i + 1
            if cf in concat:
                for j, sub in enumerate(concat[cf]):
                    values[sub] = chunk[j] if j < len(chunk) else ''
            else:
                values[cf] = chunk
            lines.append(self._render_line(card, values, fmt))
        return lines

    # ---- internals -------------------------------------------------------

    def _render_line(self, card, values, fmt):
        """Render one 80-column line from a ``{fieldname: value}`` mapping."""
        line = bytearray(b' ' * 80)
        self._place(line, 1, 6, card.ljust(6))
        for fname, spec in fmt['fields'].items():
            typestring, (start, end) = spec[0], spec[1]
            hints = _hints(spec)
            width = end - start + 1
            value = values.get(fname, '')
            if typestring in self.custom_formats:
                text = self._emit_composite(typestring, value, width)
            elif hints.get('just') == 'atomname':
                text = self._emit_atomname(value, values.get('element', ''), width)
            elif typestring == 'HxInteger':
                text = self._emit_serial(value, width)
            elif isinstance(value, list):
                text = self._emit_list(value, typestring, width)
            else:
                text = self.formatter.format(value, typestring, width, hints)
            self._place(line, start, end, text)
        return line.decode('ascii').rstrip()

    def _emit_serial(self, value, width):
        """Encode a serial-number field via the stateful decimal->hex encoder."""
        if value == '' or value is None:
            return ' ' * width
        s = self.hex(int(value))
        if len(s) > width:
            logger.warning(f'serial {int(value)} needs {len(s)} cols but the '
                           f'field is {width}; the file exceeds the hybrid-hex '
                           f'ceiling (~1M atoms) and cannot be represented')
            return s[-width:]
        return s.rjust(width)

    def _emit_list(self, items, typestring, width):
        """Join a list value into its column using the type's separator."""
        sep = LIST_SEP.get(typestring, ' ')
        s = sep.join(str(x) for x in items)
        if len(s) > width:
            logger.warning(f'list field overflows {width} cols; clipping: {s!r}')
            s = s[:width]
        return s.ljust(width)

    def _chunk(self, cf, value, fmt, concat):
        """Split a continued field's value into per-line pieces."""
        if cf in concat:                                   # e.g. SITE residues
            per = len(concat[cf])
            items = list(value) if isinstance(value, list) else []
            return [items[i:i + per] for i in range(0, len(items), per)] or [[]]
        start, end = fmt['fields'][cf][1]
        width = end - start + 1
        if isinstance(value, list):                        # WList/CList/SList
            typestring = fmt['fields'][cf][0]
            return self._pack_items(value, LIST_SEP.get(typestring, ' '), width) or [[]]
        return self._wrap_words(str(value), width) or ['']  # plain String

    @staticmethod
    def _pack_items(items, sep, width):
        """Greedily pack whole items into lines no wider than ``width``."""
        lines, cur = [], []
        for it in items:
            if cur and len(sep.join(cur + [str(it)])) > width:
                lines.append(cur)
                cur = [str(it)]
            else:
                cur.append(str(it))
        if cur:
            lines.append(cur)
        return lines

    @staticmethod
    def _wrap_words(text, width):
        """Word-wrap ``text`` at whitespace into lines no wider than ``width``."""
        lines, cur = [], ''
        for word in text.split():
            trial = (cur + ' ' + word).strip()
            if cur and len(trial) > width:
                lines.append(cur)
                cur = word
            else:
                cur = trial
        if cur:
            lines.append(cur)
        return lines

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
    'EXPDTA', 'NUMMDL', 'MDLTYP', 'AUTHOR', 'REVDAT', 'SPRSDE', 'JRNL', 'REMARK',
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


def assemble_pdb(parser, anisou=True, include_master=True, record_formats=None):
    """
    Assemble a conformant PDB document from a parsed structure.

    Emits every writable record present in ``parser.parsed`` in canonical
    section order: single-line records (types 1/3) via :meth:`PDBWriter.emit`,
    continuation/determinant-group records (types 2/4, plus ``REVDAT``) via
    :meth:`PDBWriter.emit_multiline`, and the coordinate section (``ATOM`` with
    interleaved ``ANISOU`` and chain-terminating ``TER`` cards, then
    ``HETATM``). ``REMARK`` and ``JRNL`` (type 6) are re-emitted verbatim from
    the source lines. The ``MASTER`` and ``END`` bookkeeping records are
    regenerated from the emitted content. When the input was mmCIF there are no
    source lines to pass through, so ``REMARK``/``JRNL`` are omitted and
    reported via the logger.

    Parameters
    ----------
    parser : PDBParser
        A parsed structure (``parser.parsed`` populated).
    anisou : bool, optional
        Interleave ``ANISOU`` records after their matching atom (default True).
    include_master : bool, optional
        Regenerate a ``MASTER`` record from the emitted counts (default True).
    record_formats : dict, optional
        The record-format table to write with. Defaults to the parser's active
        (dialect-aware) ``record_formats``; pass an explicit table to write a
        different dialect than the parser was built with.

    Returns
    -------
    list of str
        The document as a list of record lines (no trailing newlines).
    """
    rf = record_formats if record_formats is not None else parser.record_formats
    cf = parser.pdb_format_dict['custom_formats']
    w = PDBWriter(rf, cf)
    parsed = parser.parsed

    if 'MODEL' in parsed:
        logger.warning('multi-model coordinate sections are not reconstructed; '
                       'emitting atoms as a single flat model')

    lines = []
    counts = {}

    #: type-6 records re-emitted verbatim from the source (never re-serialized)
    passthrough_cards = {'REMARK', 'JRNL'}

    def _emit_key(key):
        if key in passthrough_cards:
            _passthrough(key)
            return
        recs = parsed.get(key)
        if recs is None:
            return
        fmt = rf.get(key, {})
        rlist = list(recs) if isinstance(recs, (list, PDBRecordList)) else [recs]
        # presence of `continues` (not the type number) is the discriminator:
        # every continuation/determinant-group record has it, including REVDAT
        # (tagged type 3 but genuinely multi-line); nothing single-line does
        if fmt.get('continues'):
            n = 0
            for r in rlist:
                out = w.emit_multiline(r, key)
                lines.extend(out)
                n += len(out)
            counts[key] = counts.get(key, 0) + n
        elif fmt.get('type') in (1, 3):
            for r in rlist:
                lines.append(w.emit(r, key))
            counts[key] = counts.get(key, 0) + len(rlist)

    def _passthrough(card):
        """Re-emit a record type's original source lines verbatim (type 6)."""
        n = 0
        for ln in getattr(parser, 'pdb_lines', []) or []:
            if ln[:6].strip() == card:
                lines.append(ln.rstrip())
                n += 1
        if n:
            counts[card] = counts.get(card, 0) + n
        elif card in {k.split('.')[0] for k in parsed}:
            logger.info(f'{card} present in parse but no source lines to pass '
                        f'through (input was not PDB); omitted from output')

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
