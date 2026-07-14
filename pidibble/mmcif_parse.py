# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
.. module:: mmcif_parse

   :synopsis: defines the MMCIF_Parser class for parsing mmCIF files

   .. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>
   
"""

from .pdbrecord import PDBRecord, PDBRecordDict, PDBRecordList
from .baserecord import BaseRecord
import logging
logger = logging.getLogger(__name__)

def split_ri(ri):
    """
    Split a residue identifier into its sequence number and insertion code.

    Parameters
    ----------
    ri : str or int
        The residue identifier, which can be a string in the format '1234A' or an integer like 1234.

    Returns
    -------
    tuple
        A tuple containing the sequence number as an integer and the insertion code as a string.
    """
    if isinstance(ri, int):  # this is no insertion code
        r = ri
        i = ''
    elif ri[-1].isdigit():  # there is no insertion code
        r = int(ri)
        i = ''
    else:
        r = int(ri[:-1])
        i = ri[-1]
    return r, i

def rectify(val):
    """
    Convert a value to its appropriate type, handling empty strings and special cases.

    Parameters
    ----------
    val : str
        The value to be rectified, which can be a string representation of a number or an empty string.

    Returns
    -------
    int or float or str
        The rectified value, which is an integer if the string represents a number, a float if it can be converted, or the original string if it cannot be converted.
    """
    if not val:
        return ''
    if val in '.?':
        return ''
    if val.isdigit():
        return int(val)
    try:
        val = float(val)
    except ValueError:
        pass
    return val

class MMCIF_Parser:
    """
    A parser for mmCIF files, handling the parsing of various formats and structures.

    Parameters
    ----------
    mmcif_formats : dict
        A dictionary defining the mmCIF formats to be parsed.
    pdb_formats : dict
        A dictionary defining the PDB formats to be parsed.
    cif_data : object
        An object containing the CIF data to be parsed.
    """
    def __init__(self, mmcif_formats, pdb_formats, cif_data):
        self.formats = mmcif_formats
        self.pdb_formats = pdb_formats
        self.global_maps = {}
        self.global_ids = {}
        self.cif_data = cif_data

    def update_maps(self, maps, row):
        """
        Update the global maps with values from a single CIF row.

        Parameters
        ----------
        maps : dict
            A dictionary of maps to update, where keys are map names and values are dictionaries with 'key' and 'value' keys.
        row : dict
            A single CIF row as an ``{attribute: value}`` dict, as returned by
            ``DataCategory.getRowAttributeDict()``.
        """
        for mapname, mapspec in maps.items():
            if not mapname in self.global_maps:
                self.global_maps[mapname] = {}
            key = rectify(row.get(mapspec['key'], ''))
            val = rectify(row.get(mapspec['value'], ''))
            if not key in self.global_maps[mapname]:
                self.global_maps[mapname][key] = val

    def update_ids(self, idmaps, row):
        """
        Update the global IDs with values from a single CIF row.

        Parameters
        ----------
        idmaps : dict
            A dictionary of ID maps, where keys are ID names and values are the corresponding CIF record field names.
        row : dict
            A single CIF row as an ``{attribute: value}`` dict, as returned by
            ``DataCategory.getRowAttributeDict()``.
        """
        for idname, idspec in idmaps.items():
            if not idname in self.global_ids:
                self.global_ids[idname] = []
            thisid = rectify(row.get(idspec, ''))
            if not thisid in self.global_ids[idname]:
                self.global_ids[idname].append(thisid)

    def _join_lookup(self, jo, on, idict):
        """
        Find the first row of category ``jo`` matching every ``on`` condition.

        Parameters
        ----------
        jo : DataCategory or None
            The category to search.
        on : list of tuple
            ``(other_attr, self_key)`` pairs; a row matches when, for each pair,
            the row's ``other_attr`` equals this record's ``self_key`` value.
        idict : dict
            The record being built, holding the already-mapped self values.

        Returns
        -------
        dict
            The matching row as an ``{attribute: value}`` dict, or ``{}`` if none.
        """
        if jo is None or not on:
            return {}
        first_attr, first_key = on[0]
        candidates = jo.selectIndices(str(idict.get(first_key, '')), first_attr)
        for idx in candidates:
            jrow = jo.getRowAttributeDict(idx)
            if all(str(jrow.get(oa, '')) == str(idict.get(sk, '')) for oa, sk in on[1:]):
                return jrow
        return {}

    def gen_dict(self, mapspec):
        """
        Generate a list of dictionaries based on the specified mapping specification.
        This method processes the mapping specification to create dictionaries that represent parsed records from the CIF data.

        Parameters
        ----------
        mapspec : dict
            A dictionary containing the mapping specification, which includes keys like 'data_obj', 'attr_map', 'splits', 'spawns_on', 'indexes', 'map_values', 'tables', 'spawn_data', 'global_maps', 'global_ids', 'list_attr', 'signal_attr', 'signal_value', 'allcaps', and 'if_dot_replace_with'.

        Returns
        -------
        list
            A list of dictionaries representing the parsed records based on the mapping specification.
        """
        idicts = []
        attr_map = mapspec.get('attr_map', {})
        splits = mapspec.get('splits', [])
        spawns_on = mapspec.get('spawns_on', None)
        indexes = mapspec.get('indexes', None)
        map_values = mapspec.get('map_values', {})
        tables = mapspec.get('tables', {})
        spawn_data = mapspec.get('spawn_data', {})
        tables = mapspec.get('tables', {})
        list_attr = mapspec.get('list_attr', {})
        sigattr = mapspec.get('signal_attr', None)
        sigval = mapspec.get('signal_value', None)
        use_signal = (sigattr is not None)
        global_maps = mapspec.get('global_maps', {})
        global_ids = mapspec.get('global_ids', {})
        spawns_on = mapspec.get('spawns_on', None)
        allcaps = mapspec.get('allcaps', [])
        if_dot_replace_with = mapspec.get('if_dot_replace_with', {})
        logger.debug(f'getting cifrec for {mapspec["data_obj"]}')
        cifrec = self.cif_data.getObj(mapspec['data_obj'])
        groupby = mapspec.get('groupby', None)
        if groupby and cifrec is not None:
            # group rows sharing an attribute value (e.g. author chain) into one
            # record and collect per-group lists — reproduces PDB grouped records
            # such as SEQRES (one record per chain, residues gathered in order).
            group_attr_map = mapspec.get('group_attr_map', {})
            collect = mapspec.get('collect', {})
            lengths = mapspec.get('lengths', {})
            groups = {}  # groupkey -> {'first': row, 'collected': {out: [vals]}}
            for idx in range(len(cifrec)):
                row = cifrec.getRowAttributeDict(idx)
                gk = row.get(groupby, '')
                if gk not in groups:
                    groups[gk] = {'first': row, 'collected': {ck: [] for ck in collect}}
                for ck, cattr in collect.items():
                    groups[gk]['collected'][ck].append(rectify(row.get(cattr, '')))
            for gk, g in groups.items():
                idict = {}
                for k, cattr in group_attr_map.items():
                    idict[k] = rectify(g['first'].get(cattr, ''))
                for ck, vals in g['collected'].items():
                    idict[ck] = vals
                for lk, ck in lengths.items():
                    idict[lk] = len(idict[ck])
                idicts.append(idict)
        elif not tables and cifrec is not None:
            # select matching rows up front (selectIndices preserves row order)
            # instead of scanning every row and testing the signal attribute;
            # signal_value may be a single value or a list of accepted values
            if use_signal:
                sigvals = sigval if isinstance(sigval, list) else [sigval]
                indices = sorted(i for sv in sigvals for i in cifrec.selectIndices(sv, sigattr))
            else:
                indices = range(len(cifrec))
            for idx in indices:
                # pull the whole row as an {attribute: value} dict once, rather
                # than issuing a positional getValue() per field
                row = cifrec.getRowAttributeDict(idx)
                if global_maps:
                    self.update_maps(global_maps, row)
                if global_ids:
                    self.update_ids(global_ids, row)
                idict = {}
                for k, v in attr_map.items():
                    if isinstance(v, dict):
                        resdict = {kk: rectify(row.get(o, '')) for kk, o in v.items()}
                        if 'resseqnumi' in resdict:
                            resdict['seqNum'], resdict['iCode'] = split_ri(resdict['resseqnumi'])
                        val = PDBRecord(resdict)
                    else:
                        val = rectify(row.get(v, ''))
                        if k == 'resseqnumi':
                            idict['seqNum'], idict['iCode'] = split_ri(val)
                        else:
                            if k in splits and ',' in val:
                                val = [rectify(x.strip()) for x in val.split(',')]
                            if k == spawns_on:
                                if isinstance(val, str) and ',' in val:
                                    val = [rectify(x.strip()) for x in val.split(',')]
                            if k in map_values:
                                mapper = self.global_maps[map_values[k]]
                                if isinstance(val, list):
                                    logger.debug(f'mapper {mapper}')
                                    logger.debug(f'list before mapping {val}')
                                    mapped_val = list(set([str(mapper[x]) for x in val]))
                                    logger.debug(f'list after mapping {mapped_val}')
                                    try:
                                        mapped_val.sort()
                                        val = mapped_val
                                    except TypeError:
                                        raise TypeError(f'could not sort list {mapped_val} at key {k}')
                                else:
                                    val = mapper[val]
                    idict[k] = val
                    if k == indexes:
                        idict['tmp_label'] = f'{k}{val}'
                for la, vn in list_attr.items():
                    from_existing = all([x in idict for x in vn])
                    if from_existing:
                        idict[la] = [idict[x] for x in vn]
                    else:
                        idict[la] = vn
                if spawns_on:
                    spdicts = self.gen_dict(mapspec['spawn_data'])
                    if isinstance(idict[spawns_on], list):
                        spawned_dicts = []
                        for v in idict[spawns_on]:
                            sd = idict.copy()
                            sd[spawns_on] = v
                            for sp in spdicts:
                                if sp['spawn_idx'] == v:
                                    break
                            else:
                                raise Exception(f'(list) cannot find spawn index for {spawns_on} = {v}; spdicts: {spdicts}')
                            spc = sp.copy()
                            del spc['spawn_idx']
                            spclabel = spc.get('tmp_label', '')
                            if 'tmp_label' in spc:
                                del spc['tmp_label']
                            sd.update(spc)
                            if 'tmp_label' in sd and spclabel != '':
                                sd['tmp_label'] = f'{sd["tmp_label"]}.{spclabel}'
                            spawned_dicts.append(sd)
                        idicts.extend(spawned_dicts)
                    else:
                        spawned_dicts = []
                        v = idict[spawns_on]
                        for sp in spdicts:
                            if sp['spawn_idx'] == v:
                                break
                        else:
                            raise Exception(f'cannot find spawn index for {spawns_on} = {v}')
                        spc = sp.copy()
                        del spc['spawn_idx']
                        spclabel = spc.get('tmp_label', '')
                        if 'tmp_label' in spc:
                            del spc['tmp_label']
                        idict.update(spc)
                        if 'tmp_label' in idict and spclabel != '':
                            idict['tmp_label'] = f'{idict["tmp_label"]}.{spclabel}'
                        idicts.append(idict)
                else:
                    idicts.append(idict)
        elif tables and cifrec is not None:
            tabledict = {}
            for tname, tspec in tables.items():
                tabledict[tname] = []
                attr_map = tspec['row_attr_map']
                bisv = tspec.get('blank_if_single_valued', [])
                for i in range(len(cifrec)):
                    row = cifrec.getRowAttributeDict(i)
                    tdict = {}
                    for k, v in attr_map.items():
                        tdict[k] = rectify(row.get(v, ''))
                        if k in bisv:
                            if len(self.global_ids[k]) < 2:
                                tdict[k] = ''
                    tabledict[tname].append(BaseRecord(tdict))
            udict = {'tables': tabledict}
            idicts.append(udict)
        # else: the mapped category is absent from this file -> emit no records

        # merge single-valued attributes drawn from other (single-row) categories,
        # e.g. CRYST1 draws cell parameters from `cell` and space group from
        # `symmetry`. Only the first row of each merged category is consulted.
        merge = mapspec.get('merge', {})
        if merge:
            for cat_name, amap in merge.items():
                mo = self.cif_data.getObj(cat_name)
                mrow = mo.getRowAttributeDict(0) if (mo is not None and len(mo)) else {}
                for idict in idicts:
                    for ak, cifattr in amap.items():
                        idict[ak] = rectify(mrow.get(cifattr, ''))

        # keyed join: for each record, look up a row of another category whose
        # attributes match this record's already-mapped values on all `match`
        # conditions ({other_attr: self_key}), and pull additional attributes
        # (scalars or nested residue dicts) from it. e.g. COMPND draws molName
        # from `entity` on {id: molID}; SHEET draws sense from struct_sheet_order
        # on {sheet_id: sheetID, range_id_2: strand}.
        join = mapspec.get('join', {})
        if join:
            for cat_name, spec in join.items():
                jo = self.cif_data.getObj(cat_name)
                on = list(spec['match'].items())  # [(other_attr, self_key), ...]
                jmap = spec['attr_map']
                for idict in idicts:
                    jrow = self._join_lookup(jo, on, idict)
                    for ak, cifattr in jmap.items():
                        if isinstance(cifattr, dict):
                            idict[ak] = PDBRecord({kk: rectify(jrow.get(o, '')) for kk, o in cifattr.items()})
                        else:
                            idict[ak] = rectify(jrow.get(cifattr, ''))

        # map literal values to replacements (e.g. SHEET sense strings
        # 'anti-parallel'/''/'parallel' -> -1/0/1 to match the PDB integer)
        value_maps = mapspec.get('value_maps', {})
        if value_maps:
            for idict in idicts:
                for k, vmap in value_maps.items():
                    if k in idict and idict[k] in vmap:
                        idict[k] = vmap[idict[k]]

        # ensure the named attributes are always lists (e.g. a single-chain
        # COMPND still yields chains=['G'] rather than 'G')
        as_list = mapspec.get('as_list', [])
        for idict in idicts:
            for k in as_list:
                v = idict.get(k, '')
                if not isinstance(v, list):
                    idict[k] = [v] if v != '' else []

        if allcaps:
            for idict in idicts:
                for k, v in idict.items():
                    if k in allcaps:
                        if isinstance(v, list):
                            idict[k] = [x.upper() if isinstance(x, str) else x for x in v]
                        elif isinstance(v, str):
                            idict[k] = v.upper()
        return idicts

    def parse(self):
        """
        Parse the mmCIF data and generate a dictionary of :class:`pdbrecord.PDBRecord` instances.
        This method processes the mmCIF formats and generates a dictionary where keys are record types and values are lists of :class:`pdbrecord.PDBRecord` instances.

        Returns
        -------
        PDBRecordDict
            A dictionary where keys are record types and values are lists of :class:`pdbrecord.PDBRecord` instances.
        """
        recdict = PDBRecordDict()
        for rectype, mapspec in self.formats.items():
            idicts = self.gen_dict(mapspec)
            for idict in idicts:
                this_key = idict.get('tmp_label', '')
                reckey = rectype if not this_key else f'{rectype}.{this_key}'
                if reckey in recdict:
                    if not isinstance(recdict[reckey], PDBRecordList):
                        recdict[reckey] = PDBRecordList([recdict[reckey]])
                    idict['key'] = reckey
                    recdict[reckey].append(PDBRecord(idict))
                else:
                    idict['key'] = reckey
                    recdict[reckey] = PDBRecord(idict)
        self._report_unmapped_categories()
        return recdict

    def _referenced_categories(self):
        """Return the set of mmCIF category names any mapspec reads from."""
        cats = set()
        for mapspec in self.formats.values():
            cats.add(mapspec.get('data_obj'))
            for directive in ('merge', 'join'):
                cats.update((mapspec.get(directive) or {}).keys())
            spawn = mapspec.get('spawn_data')
            if spawn:
                cats.add(spawn.get('data_obj'))
        cats.discard(None)
        return cats

    def _report_unmapped_categories(self):
        """
        Log a coverage summary: how many of the file's categories pidibble reads,
        and (at DEBUG) which ones it ignores. Uses the py-mmcif category index
        so users can see what data is present but not surfaced.
        """
        try:
            present = set(self.cif_data.getObjNameList())
        except Exception:
            return
        mapped = self._referenced_categories() & present
        unmapped = sorted(present - mapped)
        if unmapped:
            logger.info(f'mmCIF: read {len(mapped)} of {len(present)} categories present; '
                        f'{len(unmapped)} present but unmapped (set logging to DEBUG to list)')
            logger.debug(f'unmapped mmCIF categories: {", ".join(unmapped)}')
