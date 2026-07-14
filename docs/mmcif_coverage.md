# PDB vs. PDBx/mmCIF parsing — coverage audit & roadmap

**Status:** audit as of 2026-07-14 (pidibble 1.6.0). Developer/roadmap reference,
not end-user documentation.

## Summary

pidibble is a **near-complete PDB (legacy fixed-column) parser** and a
**deliberately minimal PDBx/mmCIF parser**. The mmCIF path extracts only enough
to feed a downstream consumer (atoms, bonds, assemblies, missing residues); it
is not a general mmCIF reader. Parsing the same entry (4TVP) both ways:

| | Record types produced |
|---|---|
| PDB | 36 |
| mmCIF | 11 (was 6 before roadmap #3) |
| in PDB but not mmCIF | 25 |

The same 4TVP mmCIF file contains **73 data categories**; pidibble maps 11 of them.

## External library

pidibble depends on **`mmcif`** (PyPI), which is the **RCSB/wwPDB `py-mmcif`**
library ("mmCIF Core Access Library", https://github.com/rcsb/py-mmcif) — the
authoritative, C++-backed reader. This is the correct choice; it is *not*
`mmcif-pdbx`/`pdbx` (a lighter third-party alternative). **No package change is
recommended.**

### API utilization: we use ~3 methods of dozens

The entire integration ([`read_mmCIF`](../pidibble/pdbparse.py) +
[`MMCIF_Parser`](../pidibble/mmcif_parse.py)) touches only:

- `dc.getObj(name)` — fetch a category
- `cat.getValue(attr, idx)` — one cell at a time, in manual positional loops
- `len(cat)` — row count

High-value unused capabilities that should shape the expansion:

| Capability | Method | Why it matters |
|---|---|---|
| Category discovery | `DataContainer.getObjNameList()` | Enumerate all categories actually present (73 in 4TVP) instead of hard-coding names; lets the nonconformance registry flag "present but unmapped." |
| Row-as-dict | `DataCategory.getRowAttributeDict(idx)`, `getRowList()`, `getAttributeList()` | Pull a whole row as `{attr: value}` in one call; collapses most of `gen_dict`'s positional index machinery and makes adding categories declarative. **Adopted 2026-07-14 (roadmap #2).** |
| Declarative selection | `selectValuesWhere(...)`, `selectIndices(...)` | Native replacement for the manual `signal_attr`/`signal_value` scan (e.g. `atom_site` by `group_PDB`, `struct_conn` by `covale`/`disulf`). **Adopted 2026-07-14 (roadmap #2).** |
| Safe optional access | `getValueOrDefault(attr, idx, default)`, `hasAttribute(name)` | Clean handling of missing/optional items instead of carefully avoiding absent attributes. |
| DDL2 dictionary | `mmcif.api.DictionaryApi` | Authoritative type/enum/mandatory info per item; could replace the heuristic `rectify()` coercion and semi-auto-generate the format map. Heavier (needs the dictionary file) — a later phase. |

**Design implication:** the 30 missing categories should be added by leaning on
`getObjNameList` + `getRowAttributeDict` + `selectValuesWhere`, not by
hand-writing 30 more positional mapspecs.

## What mmCIF currently parses

Defined in [`resources/mmcif_format.yaml`](../pidibble/resources/mmcif_format.yaml)
(7 mapspecs → 6 record keys):

| pidibble key | mmCIF category | Caveats |
|---|---|---|
| `ATOM` | `atom_site` (group_PDB=ATOM) | Full coordinate parity with PDB; the only well-tested path |
| `HETATM` | `atom_site` (group_PDB=HETATM) | Same |
| `LINK` | `struct_conn` (**covale only**) | Validated vs PDB (4TVP); `metalc`, `hydrog`, salt bridges dropped |
| `SSBOND` | `struct_conn` (disulf) | Validated vs PDB (4TVP) |
| `SEQADV` | `struct_ref_seq_dif` | Validated vs PDB (4TVP) |
| `REMARK.350` | `pdbx_struct_assembly_gen` + `pdbx_struct_oper_list` | Bio-assembly + transforms; validated vs PDB (4TVP) |
| `REMARK.465` | `pdbx_unobs_or_zero_occ_residues` | Missing residues; validated vs PDB (4TVP) |
| `HEADER` | `entry` + `struct_keywords` + `pdbx_database_status` | idCode/classification exact; `depDate` in native ISO (see below) |
| `TITLE` | `struct` | Uppercased to match PDB; validated vs PDB (4TVP) |
| `EXPDTA` | `exptl` | Experiment method; validated vs PDB (4TVP) |
| `KEYWDS` | `struct_keywords` | Uppercased list; see keywds caveat below |
| `CRYST1` | `cell` + `symmetry` | Cell/space-group/Z; validated vs PDB (4TVP) |

### Representational caveats for the new metadata records

- **`HEADER.depDate`** is exposed in mmCIF's native ISO form (`2014-06-27`),
  which differs from the PDB `DD-MON-YY` form (`27-JUN-14`). A future
  normalization pass could reconcile the two.
- **Text case:** `TITLE.title` and `KEYWDS.keywds` are uppercased so the mmCIF
  output matches the PDB convention and consumers get identical data regardless
  of source format. This intentionally discards mmCIF's richer mixed-case text.
- **`KEYWDS.keywds`** exact correspondence can be confounded by source-data
  comma placement (in 4TVP the PDB `KEYWDS` splits one phrase that the mmCIF
  `struct_keywords.text` leaves joined); this is a file-content difference, not
  a parsing bug.
- Cryo-EM entries (e.g. 8FAE) carry a dummy `cell`/`symmetry` (`a=1.0`,
  `sGroup='P 1'`, empty `z`); this is passed through as-is.

## What mmCIF does NOT parse (the 30 PDB record types)

With the mmCIF category that would supply each:

`HEADER`, `TITLE`, `EXPDTA`, `KEYWDS`, and `CRYST1` are now mapped (roadmap #3);
the remaining gaps are:

| PDB record(s) | Data | mmCIF category |
|---|---|---|
| `COMPND`, `SOURCE` | entities, chains, organism | `entity`, `entity_poly`, `entity_src_gen`, `struct_asym` |
| `REMARK 2 / 3` | resolution, refinement | `reflns`, `refine`, `refine_ls_restr` |
| `SCALE1-3` | scale matrix | `atom_sites` |
| `ORIGX1-3` | origin transform | `atom_sites` (or `database_PDB_matrix`) |
| `SEQRES` | sequence | `entity_poly_seq`, `pdbx_poly_seq_scheme` |
| `DBREF` | DB cross-refs | `struct_ref`, `struct_ref_seq` |
| `HELIX` | helices | `struct_conf` |
| `SHEET` | sheets | `struct_sheet_range`, `struct_sheet_order`, `pdbx_struct_sheet_hbond` |
| `MODRES` | modified residues | `pdbx_struct_mod_residue` |
| `SITE` | functional sites | `struct_site`, `struct_site_gen` |
| `FORMUL`, `HETNAM`, `HETSYN`, `HET` | het chemistry | `chem_comp`, `pdbx_entity_nonpoly` |
| `ANISOU` | anisotropic U | `atom_site_anisotrop` |
| `CISPEP` | cis peptides | `struct_mon_prot_cis` |
| `CONECT` | connectivity | `struct_conn`, `chem_comp_bond` |
| `REVDAT` | revision history | `pdbx_audit_revision_history` |
| `JRNL` | primary citation | `citation`, `citation_author` |
| `LINK` (non-covalent) | metal/other bonds | `struct_conn` (`metalc`, etc.) |
| `MASTER`, `END`, `TER` | bookkeeping | n/a (structural artifacts) |

## Correctness / validation status

**Resolved 2026-07-14 (roadmap #1 done).** Previously the only mmCIF↔PDB
equivalence test asserted on ATOM/HETATM only; the blocks for LINK, SSBOND,
SEQADV, REMARK.465, and REMARK.350 were commented out. They have now been
revived as five passing tests in
[`Test_mmCIF`](../tests/unit/test_rcsb.py):
`test_cif_pdb_correspondence_{links,ssbonds,seqadv,missing_residues,assembly}`.

**Why they were disabled:** the original blocks compared mmCIF **label**
numbering (`residue1` → chain A, seq 58) against PDB **author** numbering
(chain G, seq 88), which never matches. The fix is to compare mmCIF `*_auth`
residues (`residue1_auth`, `residue_auth`, `auth_*`) against PDB residues — the
same convention the atoms test already used. LINK `length` also needed a
tolerance (PDB writes 2 decimals, e.g. `1.44`, vs mmCIF `1.437`).

**Result for 4TVP:** all five categories agree element-wise with zero
mismatches — LINK 59/59, SSBOND 21/21, SEQADV 9/9, REMARK.465 69/69, and the 3
assembly transforms have matching rotation matrices and translation vectors. So
the existing mmCIF parsing of these categories is **correct**, at least for an
entry with exact PDB↔mmCIF correspondence. (Caveat: 4TVP is such an entry;
large or PDB-less entries with no legacy counterpart are not exercised by this
test.)

## Minor issues

- **Dead code** in `mmcif_parse.py`: `MMCIFDict` class and `resolve()` (a `pass`
  stub) are never instantiated or called.
- **No multi-model handling** on the mmCIF side: `pdbx_PDB_model_num` is captured
  as an atom attribute, but atoms are not grouped into `MODEL`/`ENDMDL`
  structures as in the PDB path.
- The **nonconformance registry** (1.6.0) covers only the PDB path; mmCIF coercion
  (`rectify`, `getValue`) has no equivalent reporting.

## Roadmap (priority order)

1. ~~**Revive the correspondence tests.**~~ **Done 2026-07-14.** Five tests
   revived and passing; existing mmCIF parse of LINK/SSBOND/SEQADV/REMARK.465/
   REMARK.350 verified correct against PDB for 4TVP. See "Correctness /
   validation status" above.
2. ~~**Refactor `gen_dict` onto `getRowAttributeDict` + `selectValuesWhere`.**~~
   **Done 2026-07-14.** Positional `getValue(attr, idx)` loops replaced with
   `selectIndices()` for signal filtering and `getRowAttributeDict()` for row
   access; `update_maps`/`update_ids` now take a row dict. Verified
   behavior-preserving: parsed output is byte-for-byte identical (SHA-256) for
   4tvp/8fae/4zmj, and all correspondence tests still pass. Missing attributes
   now yield `''` (via `row.get`) instead of raising `ValueError`, so mapspec
   fields absent from a given entry degrade gracefully.
3. ~~**Add header/metadata mapping.**~~ **Done 2026-07-14.** Added HEADER,
   TITLE, EXPDTA, KEYWDS, CRYST1 from `entry`/`struct`/`struct_keywords`/
   `exptl`/`cell`/`symmetry`/`pdbx_database_status`, mirroring PDB attribute
   names and validated vs PDB for 4TVP. Introduced a `merge` mapspec directive
   (pull single-valued attrs from a secondary category, e.g. CRYST1 ←
   cell+symmetry), taught `allcaps` to handle list values, and made comma
   `splits` strip whitespace. Entity/source (`COMPND`/`SOURCE` ←
   `entity`/`entity_src_gen`) deferred — its PDB token structure is a larger
   mapping job; see remaining gaps above.
4. **Add `SEQRES` equivalent** (`entity_poly_seq` / `pdbx_poly_seq_scheme`).
5. **Add secondary structure** (`struct_conf` → HELIX, `struct_sheet_range` →
   SHEET).
6. **Broaden `struct_conn`** beyond covale/disulf (at least `metalc`).
7. **Category-discovery pass** using `getObjNameList()`; have the nonconformance
   registry report categories present in the file but unmapped.
8. **Delete dead `MMCIFDict`/`resolve` code.**
9. *(Later)* Evaluate `mmcif.api.DictionaryApi` for type-correct coercion and
   format-map generation.

## References

- py-mmcif: https://github.com/rcsb/py-mmcif — API docs for `DataContainer` /
  `DataCategory`.
- PDBx/mmCIF dictionary: https://mmcif.wwpdb.org/
