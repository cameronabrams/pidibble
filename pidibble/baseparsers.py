"""

.. module:: baseparsers
   :synopsis: defines some basic string and list parsing functions
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger = logging.getLogger(__name__)

class ListParser:
    """
    A simple parser for lists of strings, with a customizable delimiter.
    """

    def __init__(self, d=','):
        self.d = d

    def parse(self, string):
        """
        Parse a string into a list of strings, using the specified delimiter.
        If no delimiter is specified, it splits on whitespace.

        Parameters
        ----------
        string : str
            The string to parse.

        Returns
        -------
        list
            A list of strings parsed from the input string.
        """
        if self.d is None:
            return [x for x in string.split() if x.strip() != '']
        else:
            return [x.strip() for x in string.split(self.d) if x.strip() != '']

def list_parse(obj, d):
    """
    A factory function to create a ListParser with a specific delimiter.

    Parameters
    ----------
    obj : type
        The class to instantiate (should be ListParser).
    d : str or None
        The delimiter to use for parsing. If None, it will split on whitespace.

    Returns
    -------
    function
        A function that takes a string and returns a list of parsed strings.
    """
    return obj(d).parse

"""
Define a dictionary of parsers for different list formats
"""
ListParsers = {
    'CList': list_parse(ListParser, ','),
    'SList': list_parse(ListParser, ';'),
    'WList': list_parse(ListParser, None),
    'DList': list_parse(ListParser, ':'),
    'LList': list_parse(ListParser, '\n')
}

_cols = """
         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890"""

class NonconformanceRegistry:
    """
    Accumulates *types* of PDB-format nonconformance encountered during a parse,
    rather than every individual instance.

    Real-world PDB files routinely deviate from the fixed-column standard (for
    example, a `charge` column that holds ``O-`` instead of a bare number). A
    single file can contain tens of thousands of instances of the same
    deviation, so logging each one is useless noise. This registry collapses
    them by *type*, keyed by ``(record_key, field, kind)`` — e.g.
    ``("ATOM", "charge", "not coercible to Float")`` — keeping a running count
    and one exemplar (the offending value and its byte range) per type so the
    end-of-parse summary can point at something concrete.
    """

    def __init__(self):
        self._types = {}  # (record_key, field, kind) -> {count, byte_range, example}

    def __len__(self):
        return len(self._types)

    def __bool__(self):
        return bool(self._types)

    def add(self, record_key, nc):
        """
        Record one nonconformance instance under its type signature.

        Parameters
        ----------
        record_key : str
            The record type the instance was found in (e.g. ``ATOM``).
        nc : dict
            A nonconformance descriptor as produced by :class:`StringParser`,
            with keys ``field``, ``kind``, ``byte_range``, and ``value``.
        """
        sig = (record_key, nc['field'], nc['kind'])
        entry = self._types.get(sig)
        if entry is None:
            self._types[sig] = {'count': 1,
                                'byte_range': nc.get('byte_range'),
                                'example': nc.get('value', '')}
        else:
            entry['count'] += 1

    def types(self):
        """
        Return the accumulated nonconformance types.

        Returns
        -------
        list of tuple
            One ``(record_key, field, kind, count, byte_range, example)`` tuple
            per distinct nonconformance type.
        """
        return [(rk, field, kind, e['count'], e['byte_range'], e['example'])
                for (rk, field, kind), e in self._types.items()]

    def report(self, log, label=''):
        """
        Emit one INFO line per nonconformance type to the given logger.

        Parameters
        ----------
        log : logging.Logger
            The logger to emit the summary to.
        label : str, optional
            A label for the source being parsed (e.g. the PDB id), included in
            the summary header.
        """
        if not self._types:
            return
        where = f' while parsing {label}' if label else ''
        log.info(f'{len(self._types)} PDB-format nonconformance type(s) encountered{where}; '
                 f'set logging to DEBUG for per-instance detail:')
        for (rk, field, kind), e in sorted(self._types.items()):
            name = f'{rk}.{field}' if field else rk
            br = e.get('byte_range')
            loc = f' at cols {br[0]}-{br[1]}' if br else ''
            eg = f' — e.g. {e["example"]!r}{loc}' if (e.get('example') or br) else ''
            log.info(f'  {name}: {e["count"]} value(s) {kind}{eg}')


class StringParser:
    """
    A parser for fixed-width strings, with a customizable field map.

    Parameters
    ----------
    fmtdict : dict
        A dictionary mapping field names to tuples of (type, byte_range).
    typemap : dict
        A dictionary mapping type names to Python types.
    allowed : dict, optional
        A dictionary mapping field values to allowed values, for validation.
    """
    def __init__(self, fmtdict, typemap, allowed={}):
        self.typemap = typemap
        self.fields = {k: v for k, v in fmtdict.items()}
        self.allowed = allowed
        # structured record of every nonconformance found by parse(); drained by
        # the caller into a NonconformanceRegistry for type-level summarizing
        self.nonconformances = []

    def parse(self, record):
        """
        Parse a fixed-width string record into a dictionary of fields.

        Parameters
        ----------
        record : str
            The fixed-width string record to parse.

        Returns
        -------
        dict
            A dictionary of fields parsed from the input record.
        """
        if len(record) > 80:
            self.nonconformances.append({'field': '', 'kind': 'record exceeds 80 bytes',
                                         'byte_range': None, 'value': record.strip()})
            logger.debug('The following record exceeds 80 bytes in length:')
            self.report_record_error(record)
            logger.debug('Stripping...')
            record = record.strip()
            if len(record) > 80:
                raise ValueError(f'Record is too long; something wrong with your PDB file?')
        input_dict = {}
        record += ' ' * (80 - len(record))  # pad
        for k, v in self.fields.items():
            # a field spec is [typestring, byte_range] with an optional third
            # element carrying writer hints (prec/just); parsing ignores it
            typestring, byte_range = v[0], v[1]
            typ = self.typemap[typestring]
            assert byte_range[1] <= len(record), f'{record} {byte_range}'
            # using columns beginning with "1" not "0"
            fieldstring = record[byte_range[0] - 1:byte_range[1]]
            fieldstring = fieldstring.rstrip()
            try:
                # if len(fieldstring)>0 and not typ==str:
                #     fieldstring=''
                input_dict[k] = '' if fieldstring == '' else typ(fieldstring)
            except (ValueError, TypeError):
                self.nonconformances.append({'field': k, 'kind': f'not coercible to {typestring}',
                                             'byte_range': byte_range, 'value': fieldstring})
                self.report_field_error(record, k)
                input_dict[k] = ''
            if typ == str:
                input_dict[k] = input_dict[k].strip()
            if fieldstring in self.allowed:
                assert input_dict[k] in self.allowed[fieldstring], f'Value {input_dict[k]} is not allowed for field {k}; allowed values are {self.allowed[fieldstring]}'
        return input_dict

    def report_record_error(self, record, byte_range=[]):
        """
        Report an error in parsing a fixed-width string record.

        Parameters
        ----------
        record : str
            The fixed-width string record that caused the error.
        byte_range : list, optional
            A list of byte ranges to highlight in the error message.
            If empty, the entire record is reported. 
        """
        if byte_range:
            record = record[:byte_range[0] - 1] + '\033[91m' + record[byte_range[0] - 1:byte_range[1]] + '\033[0m' + record[byte_range[1]:]
        repstr = _cols + '\n' + record + '|'
        logger.debug(repstr)

    def report_field_error(self, record, k):
        """
        Report an error in parsing a specific field from a fixed-width string record.

        Parameters
        ----------
        record : str
            The fixed-width string record that caused the error.
        k : str
            The field name that caused the error.
        """
        byte_range = self.fields[k][1]
        logger.debug(f'Could not parse field {k} from bytes {byte_range}:')
        self.report_record_error(record, byte_range=byte_range)

def safe_float(x):
    """
    Convert a string to a float, returning 0.0 if the string is 'nan'.
    """
    if x == 'nan':
        return 0.0
    return float(x)

def str2int_sig(arg: str):
    """
    Convert a string to an integer, returning -1 if the string is not numeric.
    If the string starts with a '-', it is returned as an integer.
    """
    stripped = arg.strip()
    if not stripped.isnumeric():
        if stripped and stripped[0] == '-':
            return int(arg)
        else:
            return -1
    return int(arg)
