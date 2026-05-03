# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

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
