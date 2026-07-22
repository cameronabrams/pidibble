# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""

.. module:: pdbparse
   :synopsis: Defines the PDBParser class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import importlib.metadata
import importlib.resources
import json
import logging
import os
import urllib.request
import yaml
from typing import Callable

import numpy as np

from mmcif.io.IoAdapterCore import IoAdapterCore
from pathlib import Path

from .baseparsers import ListParsers, ListParser, str2int_sig, safe_float, NonconformanceRegistry
from .baserecord import BaseRecordParser
from .pdbrecord import PDBRecord, PDBRecordDict, PDBRecordList
from .mmcif_parse import MMCIF_Parser
from .hex import AtomSerialParser, str2atomSerial, hex_reset

logger = logging.getLogger(__name__)
__version__ = importlib.metadata.version("pidibble")

class PDBParser:
    """
    A class for parsing PDB files and extracting structured data.
    This class handles fetching PDB files, reading them, and parsing their contents into structured records
    based on predefined formats.

    Attributes
    ----------

    parsed : PDBRecordDict
        A dictionary containing parsed records, where keys are record types and values are :class:`.pdbrecord.PDBRecord` instances or lists of instances.
        This dictionary is populated after parsing the PDB or mmCIF file.

    mappers : dict
        A dictionary of mappers for parsing different data types, including custom formats and delimiters.

    pdb_lines : list
        A list of lines read from the PDB file. Empty if no input file is provided.

    cif_data : dict
        A dictionary containing the parsed mmCIF data. Empty if no input file is provided.
    """
    def __init__(self,
                 input_format: str = 'PDB',
                 overwrite: bool = False,
                 source_db: str = None,
                 source_id: str = None,
                 filepath: str | Path = None,
                 mappers: dict[str, Callable] = None,
                 comment_chars: list[str] = ['#'],
                 dialect: str = 'standard',
                 pdb_format_file: str = 'pdb_format.yaml',
                 mmcif_format_file: str = 'mmcif_format.yaml',
                 **kwargs):
        logger.debug(f'Pidibble v. {__version__}')
        self.input_format = input_format
        self.overwrite = overwrite
        if 'PDBcode' in kwargs:
            source_db = source_db or 'rcsb'
            source_id = source_id or kwargs['PDBcode']
        self.source_db = source_db
        self.source_id = source_id
        self.filepath = Path(filepath) if filepath else None
        self._atom_serial_parser = AtomSerialParser()
        self.mappers = {
            'HxInteger': self._atom_serial_parser,
            'Integer': str2int_sig,
            'String': str,
            'Float': safe_float,
        }
        if mappers is not None:
            self.mappers.update(mappers)
        self.mappers.update(ListParsers)
        self.comment_chars = comment_chars
        self.dialect = dialect
        self.pdb_lines = []
        self.cif_data = {}

        self.parsed = PDBRecordDict()
        self.nonconformances = NonconformanceRegistry()
        self.pdb_format_dict = self._load_format(pdb_format_file)
        self.mmcif_format_dict = self._load_format(mmcif_format_file)
        # the active record-format table for this instance's dialect; both
        # parsing and writing read this so they stay exact inverses
        self.record_formats = self._apply_dialect(self.pdb_format_dict, dialect)

        # update mappers with delimiters and custom formats
        delimiter_dict = self.pdb_format_dict.get('delimiters', {})
        for map, d in delimiter_dict.items():
            if not map in self.mappers:
                self.mappers[map] = ListParser(d).parse
        cformat_dict = self.pdb_format_dict.get('custom_formats', {})
        for cname, cformat in cformat_dict.items():
            if not cname in self.mappers:
                self.mappers[cname] = BaseRecordParser(cformat, self.mappers).parse

    @staticmethod
    def _apply_dialect(pdb_format_dict, dialect):
        """
        Build the record-format table for a dialect.

        Returns a copy of the base ``record_formats`` with the dialect's
        coordinate-record overrides merged in. ``'standard'`` is strict wwPDB;
        ``'charmm'`` swaps ``ATOM``/``HETATM``/``TER`` for the wide-resName,
        segID-bearing CHARMM layouts (``charmm_formats`` in the YAML), keeping
        x/y/z pinned at columns 31-54.

        Parameters
        ----------
        pdb_format_dict : dict
            The loaded PDB format file.
        dialect : str
            ``'standard'`` (default) or ``'charmm'``.

        Returns
        -------
        dict
            The active ``{record_key: format}`` mapping for the dialect.
        """
        record_formats = dict(pdb_format_dict['record_formats'])
        if dialect == 'charmm':
            record_formats.update(pdb_format_dict.get('charmm_formats', {}))
        elif dialect != 'standard':
            raise ValueError(f"unknown dialect {dialect!r}; expected 'standard' or 'charmm'")
        return record_formats

    @staticmethod
    def _load_format(filename: str) -> dict:
        """Load a YAML format file, checking CWD first then the package resources."""
        local = Path(filename)
        if local.is_file():
            logger.debug(f'Pidibble uses local config file {local}')
            return yaml.safe_load(local.read_text())
        resource = importlib.resources.files('pidibble') / 'resources' / filename
        try:
            text = resource.read_text(encoding='utf-8')
        except FileNotFoundError:
            raise FileNotFoundError(
                f'{filename} not found in CWD ({Path.cwd()}) or package resources'
            )
        logger.debug(f'Pidibble uses installed config file {filename}')
        return yaml.safe_load(text)

    def fetch(self):
        """
        Fetch the PDB file based on the provided PDB code or AlphaFold ID.
        This method checks if the PDB code or AlphaFold ID is provided, constructs the appropriate file path,
        and attempts to download the file from the PDB or AlphaFold API.

        Returns
        -------
        bool
            True if the file was successfully fetched, False otherwise.
        """
        assert self.source_db is not None or self.filepath is not None, f'source_db {self.source_db} or filepath {self.filepath} must be specified for fetch()'
        if self.source_db is not None and self.source_id is None:
            raise ValueError(f'You must specify a source ID code for source_db {self.source_db}')
        if self.source_db is not None and self.source_db not in ['rcsb', 'alphafold', 'opm']:
            raise ValueError(f'Source db {self.source_db} is not recognized.')

        if self.filepath is not None:
            if not self.filepath.exists():
                raise FileNotFoundError(f'{self.filepath.name} not found.')
            return True

        match self.source_db:
            case 'rcsb':
                if self.input_format == 'PDB':
                    self.filepath = Path(f'{self.source_id}.pdb')
                elif self.input_format == 'mmCIF':
                    self.filepath = Path(f'{self.source_id}.cif')
                else:
                    logger.warning(f'Input format {self.input_format} not recognized; using PDB')
                    self.filepath = Path(f'{self.source_id}.pdb')
                BASE_URL = self.pdb_format_dict['BASE_URL']
                target_url = os.path.join(BASE_URL, self.filepath.name)
                if not self.filepath.exists() or self.overwrite:
                    try:
                        urllib.request.urlretrieve(target_url, self.filepath.name)
                    except (urllib.error.URLError, OSError):
                        logger.warning(f'Could not fetch {self.filepath.name} from {self.source_db}')
                        return False
                return True
            case 'alphafold':
                self.filepath = Path(f'{self.source_id}.pdb')
                BASE_URL = self.pdb_format_dict['ALPHAFOLD_API_URL']
                target_url = os.path.join(BASE_URL, self.source_id)
                try:
                    urllib.request.urlretrieve(target_url + r'?key=' + self.pdb_format_dict['ALPHAFOLD_API_KEY'], f'{self.source_id}.json')
                except (urllib.error.URLError, OSError):
                    logger.warning(f'Could not fetch metadata for entry with accession code {self.source_id} from AlphaFold')
                    return False
                with open(f'{self.source_id}.json') as f:
                    result = json.load(f)
                try:
                    urllib.request.urlretrieve(result[0]['pdbUrl'], self.filepath.name)
                except (urllib.error.URLError, OSError):
                    logger.warning(f'Could not retrieve {result[0]["pdbUrl"]}')
                    return False
                return True
            case 'opm':
                self.filepath = Path(f'{self.source_id}.pdb')
                BASE_URL = self.pdb_format_dict['OPM_URL']
                target_url = os.path.join(BASE_URL, self.filepath.name)
                if not self.filepath.exists() or self.overwrite:
                    try:
                        urllib.request.urlretrieve(target_url, self.filepath.name)
                    except (urllib.error.URLError, OSError):
                        logger.warning(f'Could not fetch {self.filepath.name} from {self.source_db}')
                        return False
                    logger.warning(f'Stripping blanks and END lines from OPM pdb')
                    badlines = self.filepath.read_text().split('\n')
                    with open(self.filepath.name, 'w') as f_base:
                        with open(f'{self.filepath.stem}-dum.pdb', 'w') as f_dum:
                            f = f_base
                            for line in badlines:
                                sline = line.strip()
                                if not sline.startswith('END') and len(sline) > 0:
                                    f.write(sline + '\n')
                                if sline.startswith('END') and f is f_base:
                                    f = f_dum
                    logger.debug(f'Generated {self.filepath.name} and {self.filepath.stem}-dum.pdb')

                return True
            case '_':
                logger.debug(f'Source db {self.source_db} is not recognized.')
                return False

    def read_PDB(self):
        """
        Read the PDB file and store its lines in :attr:`PDBParser.pdb_lines`.
        This method opens the PDB file, reads its contents, and splits it into lines.
        If the last line is empty, it removes it from the list of lines.
        """
        with open(self.filepath, 'r') as f:
            self.pdb_lines = f.read().split('\n')
            if self.pdb_lines[-1] == '':
                self.pdb_lines = self.pdb_lines[:-1]

    def read_mmCIF(self):
        """
        Read the mmCIF file and store its data in :attr:`PDBParser.cif_data`.
        This method uses the :class:`mmcif.io.IoAdapterCore.IoAdapterCore` to read the
        mmCIF file and store the data in :attr:`PDBParser.cif_data`.
        """
        io = IoAdapterCore()
        l_dc = io.readFile(self.filepath)
        self.cif_data = l_dc[0]

    def read(self):
        """
        Read the PDB or mmCIF file based on the input format.
        This method checks the input format and calls the appropriate read method.
        """
        if self.input_format == 'mmCIF':
            self.read_mmCIF()
        else:
            self.read_PDB()

    def parse_base(self):
        """
        Parse the base records from the PDB or mmCIF file.
        This method initializes the parsing process based on the input format.
        If the input format is mmCIF, it uses the :class:`.mmcif_parse.MMCIF_Parser` to parse the mmCIF data.
        If the input format is PDB, it uses the :class:`.pdbrecord.PDBRecord` class to parse the PDB lines.
        """
        if self.input_format == 'mmCIF':
            self.parse_mmCIF()
        else:
            self.parse_PDB()

    def parse_mmCIF(self):
        """
        Parse the mmCIF data and generate a dictionary of :class:`.pdbrecord.PDBRecord` instances.
        This method uses the :class:`.mmcif_parse.MMCIF_Parser` to parse the mmCIF data and store the parsed records
        in :attr:`PDBParser.parsed`.
        """
        mmcif_parser = MMCIF_Parser(self.mmcif_format_dict, self.record_formats, self.cif_data)
        self.parsed = mmcif_parser.parse()

    def parse_PDB(self):
        """
        Parse the PDB lines and generate a dictionary of :class:`.pdbrecord.PDBRecord` instances.
        This method iterates through the PDB lines, identifies the record type based on the first character,
        and creates a new :class:`.pdbrecord.PDBRecord` instance for each record.
        It handles different record types, including continuation records and grouped records.
        """
        self._atom_serial_parser.reset()
        self.nonconformances = NonconformanceRegistry()
        record_formats = self.record_formats
        key = ''
        record_format = {}
        group_open_record = None
        for i, pdbrecord_line in enumerate(self.pdb_lines):
            tc = pdbrecord_line[0]
            if tc in self.comment_chars:
                continue
            pdbrecord_line += ' ' * (80 - len(pdbrecord_line))
            base_key = pdbrecord_line[:6].strip()
            assert base_key in record_formats, f'{base_key} is not found in among the available record formats'
            base_record_format = record_formats[base_key]
            record_type = base_record_format['type']
            new_record = PDBRecord.newrecord(base_key, pdbrecord_line, base_record_format, self.mappers, registry=self.nonconformances)
            key = new_record.key
            record_format = new_record.format
            if record_type in [1, 2, 6]:
                if not key in self.parsed:
                    self.parsed[key] = new_record
                else:
                    # this must be a continuation record
                    assert record_type != 1, f'{key} may not have continuation records'
                    root_record = self.parsed[key]
                    root_record.continue_record(new_record, record_format, all_fields=('REMARK' in key))
            elif record_type in [3, 4, 5]:
                if not key in self.parsed:
                    # this is necessarily the first occurance of a record with this key, but since there can be multiple instances this must be a list of records
                    if 'groupuntil' in record_format:
                        group_open_record = new_record
                        logger.debug(f'opening group {group_open_record.serial} until {group_open_record.format["groupuntil"]}')
                    if group_open_record is not None and key == group_open_record.format['groupuntil']:
                        logger.debug(f'closing group {group_open_record.serial}')
                        group_open_record = None
                    if 'groupby' in record_format:
                        tok = new_record.format['groupby'].split('.')
                        if group_open_record is not None:
                            if tok[0] == group_open_record.key:
                                groupid = getattr(group_open_record, tok[1])
                                setattr(new_record, group_open_record.key.lower(), groupid)
                    self.parsed[key] = PDBRecordList([new_record])
                else:
                    # this is either
                    # (a) a continuation record of a given key.(determinants)
                    # or
                    # (b) a new set of (determinants) on this key
                    # note (b) is only option if there are no determinants
                    # first, look for key.(determinants)
                    root_record = None
                    if 'determinants' in record_format:
                        nrd = [new_record.__dict__[k] for k in record_format['determinants']]
                        for r in self.parsed[key]:
                            td = [r.__dict__[k] for k in record_format['determinants']]
                            if nrd == td:
                                root_record = r
                                break
                    if root_record:
                        # case (a)
                        assert root_record.continuation < new_record.continuation, f'continuation parsing error {record_type}'
                        root_record.continue_record(new_record, record_format)
                    else:
                        # case (b)
                        if 'groupuntil' in record_format:
                            group_open_record = new_record
                            logger.debug(f'opening group {group_open_record.serial} until {group_open_record.format["groupuntil"]}')
                        if group_open_record is not None and key == group_open_record.format['groupuntil']:
                            logger.debug(f'closing group {group_open_record.serial}')
                            group_open_record = None
                        if 'groupby' in record_format:
                            tok = new_record.format['groupby'].split('.')
                            if group_open_record is not None:
                                if tok[0] == group_open_record.key:
                                    groupid = getattr(group_open_record, tok[1])
                                    setattr(new_record, group_open_record.key.lower(), groupid)
                        self.parsed[key].append(new_record)

    def post_process(self):
        """
        Post-process the parsed records to handle embedded records, tokens, and tables.
        This method checks if the input format is mmCIF and processes the records accordingly.
        If the input format is PDB, it processes the records to handle embedded records, tokens, and tables.
        """
        if self.input_format != 'mmCIF':
            self.parse_embedded_records()
            self.parse_tokens()
            self.parse_tables()

    # def parse_models(self):
    #     n_models=self.parsed.get('NUMMDL',1)
    #     for i in range(n_models):
    #         self.parsed['MODEL'][i+1]={}
    #         # in progress

    def parse_embedded_records(self):
        """
        Parse embedded records within the parsed records.
        This method iterates through the parsed records and checks if any record has embedded records.
        If an embedded record is found, it calls the :meth:`.pdbrecord.PDBRecord.parse_embedded` method to parse the embedded records.
        It updates the :attr:`PDBParser.parsed` dictionary with the new parsed records.
        """
        new_parsed_records = {}
        for key, p in self.parsed.items():
            if isinstance(p, PDBRecord):
                rf = p.format
                if 'embedded_records' in rf:
                    new_parsed_records.update(p.parse_embedded(self.record_formats, self.mappers))
            elif isinstance(p, PDBRecordList):
                for q in p:
                    rf = q.format
                    if 'embedded_records' in rf:
                        new_parsed_records.update(q.parse_embedded(self.record_formats, self.mappers))
        self.parsed.update(new_parsed_records)

    def parse_tokens(self):
        """
        Parse tokens within the parsed records.
        This method iterates through the parsed records and checks if any record has token formats.
        If a token format is found, it calls the :meth:`.pdbrecord.PDBRecord.parse_tokens` method to parse the tokens.
        It updates the :attr:`PDBParser.parsed` dictionary with the new parsed records.
        """
        for key, p in self.parsed.items():
            if isinstance(p, PDBRecord):
                rf = p.format
                if 'token_formats' in rf:
                    p.parse_tokens(self.mappers)
            elif isinstance(p, PDBRecordList):
                for q in p:
                    rf = q.format
                    if 'token_formats' in rf:
                        q.parse_tokens(self.mappers)

    def parse_tables(self):
        """
        Parse tables within the parsed records.
        This method iterates through the parsed records and checks if any record has table formats.
        If a table format is found, it calls the :meth:`.pdbrecord.PDBRecord.parse_tables` method to parse the tables.
        It updates the :attr:`PDBParser.parsed` dictionary with the new parsed records.
        """
        for key, p in self.parsed.items():
            if isinstance(p, PDBRecordList):
                continue  # don't expect to read a table from a multiple-record entry
            rf = p.format
            if 'tables' in rf:
                p.parse_tables(self.mappers)

    def parse(self):
        """
        Parse the PDB or mmCIF file and generate a dictionary of :class:`.pdbrecord.PDBRecord` instances.
        This method first fetches the PDB or mmCIF file based on the provided PDB code or AlphaFold ID.
        It then reads the file and parses its contents into structured records.
        If the input format is mmCIF, it uses the :class:`.mmcif_parse.MMCIF_Parser` to parse the mmCIF data.
        If the input format is PDB, it uses the :class:`.pdbrecord.PDBRecord` class to parse the PDB lines.

        Returns
        -------
        self : PDBParser
            The instance of :class:`.pdbrecord.PDBRecord` containing the parsed records.
        """
        if self.fetch():
            self.read()
            self.parse_base()
            self.post_process()
            self.nonconformances.report(logger, label=self.source_id or (self.filepath.name if self.filepath else ''))
        else:
            logger.warning(f'No data.')
        return self

    def write_PDB(self, filename: str = None, anisou: bool = True, include_master: bool = True,
                  dialect: str = None):
        """
        Write the parsed structure back out as a conformant PDB file.

        This assembles every writable record in canonical section order —
        single-line records (types 1/3), continuation and determinant-group
        records (types 2/4, e.g. ``TITLE``, ``COMPND``, ``SEQRES``, ``SITE``),
        and ``TER`` — reconstructs the coordinate section, passes ``REMARK`` and
        ``JRNL`` through verbatim from the source, and regenerates the
        ``MASTER``/``END`` bookkeeping records. When the parse came from mmCIF
        (no PDB source lines), ``REMARK``/``JRNL`` are omitted and reported to
        the logger.

        Parameters
        ----------
        filename : str, optional
            Destination path. If omitted, no file is written.
        anisou : bool, optional
            Interleave ``ANISOU`` records after their atoms (default True).
        include_master : bool, optional
            Regenerate a ``MASTER`` record from the emitted counts (default True).
        dialect : str, optional
            Column dialect to write, ``'standard'`` or ``'charmm'``. Defaults to
            the dialect this parser was constructed with. The CHARMM dialect
            widens ``resName`` to 6 columns and writes the authoritative segID
            column (73-76) while pinning x/y/z at columns 31-54.

        Returns
        -------
        list of str
            The assembled document as a list of record lines.
        """
        from .pdbwrite import assemble_pdb
        record_formats = self._apply_dialect(self.pdb_format_dict, dialect) if dialect else None
        lines = assemble_pdb(self, anisou=anisou, include_master=include_master,
                             record_formats=record_formats)
        if filename:
            with open(filename, 'w') as f:
                f.write('\n'.join(lines) + '\n')
        return lines

def get_symm_ops(rec: PDBRecord):
    """
    Extract the symmetry operations from a PDB record.
    This function processes the symmetry operations from a PDB record and returns the transformation matrix and translation vector.

    Parameters
    ----------
    rec : :class:`.pdbrecord.PDBRecord`
        The PDBRecord instance containing the symmetry operations.

    Returns
    -------
    M : :class:`numpy.ndarray`
        The 3x3 transformation matrix.
    T : :class:`numpy.ndarray`
        The 3x1 translation vector.
    """
    M = np.identity(3)
    T = np.array([0., 0., 0.])
    if not (hasattr(rec, 'row') and hasattr(rec, 'coordinate')):
        raise ValueError('Invalid PDBRecord: missing row or coordinate attributes')
    assert len(rec.row) == 3, f'a transformation matrix record should not have more than 3 rows'
    for c, r in zip(rec.coordinate, rec.row):
        row = c - 1
        M[row][0] = r.m1
        M[row][1] = r.m2
        M[row][2] = r.m3
        T[row] = r.t
    return M, T
