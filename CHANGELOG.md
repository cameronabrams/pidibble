# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- PDB *writing*: parsed structures can be serialized back to conformant
  fixed-column PDB. `PDBParser.write_PDB()` assembles a document in canonical
  section order, reconstructs the coordinate section (`ATOM` with interleaved
  `ANISOU` and chain-terminating `TER` cards, then `HETATM`), and regenerates
  the `MASTER`/`END` bookkeeping records from the emitted content. The
  record-level engine (`pidibble.pdbwrite.PDBWriter`) is the inverse of the
  parser, driven by the same field specs plus optional per-field writer hints
  (`{prec, just}`) carried as a third element in the YAML field definitions.
- Coverage spans all four writable record families: single-line records
  (types 1/3, plus `TER`), continuation records (type 2 — `TITLE`, `COMPND`,
  `SOURCE`, `KEYWDS`, `AUTHOR`, …), and determinant-group records (type 4 —
  `SEQRES`, `HETNAM`, `HETSYN`, `FORMUL`, `SITE`, and the multi-line `REVDAT`),
  re-wrapped/chunked across numbered continuation lines. `REMARK` and `JRNL`
  (type 6) are re-emitted verbatim from the source lines.
- A full `parse -> write -> re-parse` round-trip preserves every parsed record
  type and all field values on 4ZMJ (60 keys) and 4TVP (64 keys, incl. `SITE`);
  the regenerated `MASTER` matches the original entry's byte-for-byte, and the
  coordinate/`SEQRES`/`HETNAM`/`FORMUL`/`KEYWDS`/`TITLE` records re-serialize
  byte-exactly.
- Hexadecimal serial numbers for structures with more than 99999 atoms (e.g.
  large solvated systems): `HexSerialEncoder` is the exact inverse of the
  parser's `AtomSerialParser`, switching to hex once a serial passes 99999 and
  staying hex thereafter — including small `CONECT` back-references, which the
  parser reads as hex once tripped. Serials round-trip up to the 5-column
  hybrid-hex ceiling (`0xFFFFF` = 1 048 575 atoms).

### Changed
- Field specs may now carry an optional third element with writer formatting
  hints; parsing ignores it (the byte-range unpacking is now width-tolerant), so
  the change is fully backward-compatible.

### Not yet supported
- Re-serialization of `REMARK`/`JRNL` from the parsed model (they are passed
  through from the source instead, so they are omitted when the input was
  mmCIF), and multi-model coordinate sections.

## [1.7.2] - 2026-07-21

### Added
- Continuous-integration workflow (`tests.yaml`) running the unit tests on
  Python 3.10–3.12 and the documentation doctests on every push and pull
  request. Previously tests ran only on tagged releases.
- README now documents PDBx/mmCIF parsing (with a worked example) and carries
  version, Python-versions, license, tests, docs, and downloads badges.

### Changed
- Documentation overhauled from a single quickstart page into a multi-page User
  Guide — loading structures, the parsed data model, working with records,
  biological assemblies and symmetry, PDBx/mmCIF, large structures, nonconformant
  files, and advanced/customization — plus a supported-record-types coverage
  table. Numerous example corrections; every doctest is now validated
  (`sphinx -b doctest`), with `NORMALIZE_WHITESPACE` enabled so pretty-printed
  output validates on content rather than incidental spacing. The docs landing
  page is retitled from "Welcome to Pidibble's documentation!" to "pidibble".

## [1.7.1] - 2026-07-15

### Added
- mmCIF `REMARK.350` (biological assembly) records now expose `header_label`, the raw label asym_id list, alongside the existing `header` (author chains, mapped via the chainmap). This preserves the full per-chain (label) detail — including glycan/ligand chains — that the deduplicated author-chain list collapses.

## [1.7.0] - 2026-07-14

### Added
- Substantially expanded PDBx/mmCIF parsing: coverage grew from 6 to 16 record types (reading 22 of a typical file's ~73 categories). New mmCIF-mapped records: `HEADER`, `TITLE`, `EXPDTA`, `KEYWDS`, `CRYST1`, `SEQRES`, `HELIX`, `SHEET` (complete — ranges, sense, numStrands, H-bond registration), and `COMPND`/`SOURCE` (as native per-entity records). Metal coordination (`metalc`) now folds into `LINK`. Every mapping is validated against the legacy PDB parse.
- mmCIF↔PDB correspondence test suite revived and extended (atoms, links, ssbonds, seqadv, missing residues, assembly, metadata, seqres, helix, sheet, entities, metal coordination).
- mmCIF coverage summary logged at parse time: how many of a file's categories pidibble reads, with the unmapped ones listed at DEBUG level.
- `docs/mmcif_coverage.md`: PDB vs PDBx/mmCIF coverage audit and roadmap.
- Reusable declarative mmCIF mapspec directives: `merge`, `join` (multi-key `match:` conditions with nested residue dicts), `groupby`/`collect`/`lengths`, `as_list`, `value_maps`, list-aware `allcaps`, and multi-value `signal_value`.

### Changed
- mmCIF `gen_dict` refactored onto the py-mmcif higher-level API (`selectIndices`, `getRowAttributeDict`) instead of positional `getValue` loops; verified byte-for-byte identical output on existing records.
- A mapped-but-absent mmCIF category now emits no records (previously produced a spurious empty `{'tables': {}}` record).

### Removed
- Dead `MMCIFDict` class and `resolve()` stub in `mmcif_parse.py`.

## [1.6.0] - 2026-07-14

### Added
- `NonconformanceRegistry`: PDB-format deviations (e.g. a `charge` column that cannot be coerced to its declared type, or a record exceeding 80 bytes) are now aggregated by *type* — keyed by `(record_key, field, kind)` — rather than logged per instance. At the end of a parse, `PDBParser` emits a single INFO summary line per type with a count and one exemplar, replacing what could be tens of thousands of repeated warnings. The registry is available programmatically as `PDBParser.nonconformances` (`.types()`, `len()`, `bool()`).

### Changed
- Per-instance field/record parse-error messages demoted from WARNING to DEBUG; the INFO type-summary is the default-visible output, with full per-instance detail available by setting the logger to DEBUG.

## [1.5.4] - 2026-07-14

### Fixed
- ATOM/HETATM/ANISOU `charge` field retyped from `Float` to `String`. The PDB charge column (79–80) is a formatted string (e.g. `1-`, `2+`), never a plain float, so every populated charge previously failed `float()` and emitted a spurious "Could not parse field charge" warning while discarding the value. Charges are now captured correctly (e.g. `1-`).
- Off-by-one in `StringParser.report_record_error()` highlight: the error display dropped the character at the field's first column and coloured a window shifted one column to the right of the offending field. The highlight now wraps exactly the columns given by the field's byte range.

## [1.5.3] - 2026-05-03

### Fixed
- `PDBParser(PDBcode=...)` legacy parameter silently swallowed by `**kwargs`; restored via explicit mapping to `source_id`/`source_db='rcsb'`
- Eight bare `except:` clauses narrowed to specific exception types throughout
- Bitwise `&=` used as logical AND in `BaseRecord.empty()` replaced with `all()`
- `str2int_sig` raised `IndexError` on empty-string input; added guard before index access
- `!= None` / `== None` comparisons replaced with `is not None` / `is None` throughout
- `type(x) == T` identity comparisons replaced with `isinstance(x, T)` throughout
- `elif type(p) == list:` in `parse_tokens()` was unreachable; corrected to `isinstance(p, PDBRecordList)`
- Mutable default arguments `hold={}` and `hold=[]` in `gather_token` and `header_check` replaced with `None`

### Changed
- Global hex-tripped flag in `hex.py` replaced with `AtomSerialParser` callable class; each `PDBParser` instance owns its own state, eliminating thread-safety hazard and cross-parse contamination
- `mappers` default argument changed from a mutable dict literal to `None`; default built inside `__init__`
- Resource YAML files now loaded via `importlib.resources.files()` instead of `os.path.dirname(resources.__file__)`, ensuring correct behavior inside wheels and zip archives
- `resources/__init__.py` removed; `resources/` is now a plain data directory, not a sub-package
- `MANIFEST.in` removed; hatchling auto-discovers package data from git-tracked files
- PEP 8 spacing applied throughout all source files
- `parse_embedded()` refactored: `triggered`/`capturing` boolean pair replaced with explicit `_EmbedState` enum (`SEARCHING`, `PRE_CAPTURE`, `CAPTURING`); setup phase extracted into `_setup_embed_context()`

### Added
- `CHANGELOG.md` with full release history
- `scripts/release.sh` for automated version rotation, tagging, and push
- `[test]` optional dependency group in `pyproject.toml`

---

## [1.5.2] - 2025-09-16

### Added
- Capability to download structures from OPM and split DUM residues into a separate PDB file

## [1.5.1] - 2025-09-15

### Fixed
- Minor fixes following 1.5.0

## [1.5.0] - 2025-09-15

### Added
- Initial OPM support

## [1.4.2] - 2025-08-11

### Changed
- Updated AlphaFold interface to current API

## [1.4.1] - 2025-08-07

### Fixed
- Parsing bug in `PDBRecordList`

## [1.4.0] - 2025-07-29

### Added
- `PDBRecordList` and `PDBRecordDict` classes

## [1.3.3] - 2025-07-29

### Fixed
- Bugs where missing records were incorrectly assumed present during mmCIF parsing

## [1.3.2] - 2025-07-25

### Added
- `filepath` parameter in `PDBParser()` for transparent reading of local files

## [1.3.1] - 2025-07-25

### Fixed
- Minor fixes following 1.3.0

## [1.3.0] - 2025-07-16

### Changed
- Streamlined class attribute usage throughout
- Full API documentation published

## [1.2.3] - 2025-03-06

### Fixed
- Negative residue sequence numbers now parsed correctly

## [1.2.2] - 2025-03-04

### Fixed
- Minor fixes following 1.2.1

## [1.2.1] - 2024-10-01

### Fixed
- Hexadecimal serial number parsing issues (again)

## [1.2.0] - 2024-09-08

### Fixed
- Hexadecimal serial number parsing issues

## [1.1.9] - 2024-08-08

### Fixed
- `nan` values and `*` filler characters in numeric fields now handled gracefully

## [1.1.8] - 2024-07-15

### Added
- Unstructured REMARK records (e.g. from PACKMOL output) parsed as `REMARK.-1`

## [1.1.6] - 2024-07-11

### Fixed
- Hexadecimal atom serial detection when no `a-f` characters are present, based on value exceeding 99999

## [1.1.5] - 2024-07-11

### Fixed
- Hex-or-integer detection now restricted to atom serial number fields only

## [1.1.4] - 2024-07-11

### Added
- Hexadecimal atom serial number support for files with more than 99999 atoms

## [1.1.3] - 2024-03-21

### Added
- Ability to group records into models for multi-model PDB entries

## [1.1.2] - 2024-02-28

### Added
- Ability to fetch structures from the AlphaFold database

## [1.1.1] - 2023-09-19

### Added
- Version detection via `importlib.metadata`

## [1.0.9.1] - 2023-08-28

### Added
- Limited mmCIF parsing: `ATOM`, `HETATM`, `SSBOND`, `LINK`, `SEQADV`, `REMARK 350`, and `REMARK 465` records

## [1.0.8] - 2023-08-23

### Fixed
- Variations in how symmetry operation matrices are represented in PDB files

## [1.0.7.7] - 2023-08-22

### Changed
- Cleaned up logging throughout

## [1.0.7.6] - 2023-08-22

### Fixed
- Leading whitespace in `resName` field of `Residue10` record sometimes ignored

## [1.0.7.5] - 2023-08-18

### Added
- Support for four-letter residue names

## [1.0.7.4] - 2023-08-04

### Added
- Logging functionality

## [1.0.7.3] - 2023-08-03

### Changed
- Improved parsing of BIOMT transforms

## [1.0.7.2] - 2023-08-03

### Added
- Documentation stub on ReadTheDocs

## [1.0.7.1] - 2023-08-03

### Added
- Support for split BIOMT tables and REMARK 280, 375, 650, and 700

## [1.0.7] - 2023-08-01

### Added
- Pretty-print support for parsed records

## [1.0] - 2023-07-01

### Added
- Initial release
