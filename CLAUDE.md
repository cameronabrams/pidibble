# CLAUDE.md — pidibble

Working notes for anyone (human or agent) maintaining this repo. Public and
tracked; keep machine-specific paths out of it.

## What pidibble is

A parser — and, since 1.8.0, a writer — for Protein Data Bank structures in the
legacy fixed-column PDB format (v3.3, 2011) and in PDBx/mmCIF. Unlike most PDB
readers it yields meaningfully parsed objects for *all* standard record types,
not just `ATOM`/`CONECT`. The mmCIF path is deliberately narrower: it maps a
subset of categories onto the *same* record objects so downstream code does not
branch on file format (`docs/mmcif_coverage.md` is the coverage audit).

Downstream consumer: **pestifer**, which delegates structure I/O to pidibble and
pins a floor version. Coordinate before releasing — a bump can move that floor.

## Layout

| Path | What |
|---|---|
| `pidibble/pdbparse.py` | `PDBParser` — fetch, read, parse, `citations()`, `write_PDB()` |
| `pidibble/pdbrecord.py` | `PDBRecord`, embedded/token/table sub-parsing |
| `pidibble/baserecord.py`, `baseparsers.py` | field parsers, `NonconformanceRegistry` |
| `pidibble/mmcif_parse.py` | mmCIF → the same record objects |
| `pidibble/pdbwrite.py` | `PDBWriter`, `assemble_pdb` — the inverse of the parser |
| `pidibble/citation.py` | `Citation`, `CitationId`, RCSB enrichment |
| `pidibble/hex.py` | decimal→hex atom serials (`AtomSerialParser`, `HexSerialEncoder`) |
| `pidibble/resources/*.yaml` | **the format specs — see below** |
| `docs/source/` | Sphinx docs; `guide/` is user-facing, `api/` is autodoc |
| `docs/mmcif_coverage.md`, `docs/testing-notes.md` | developer references, dated |

## Tests

The test path is **`tests/unit`**. `tests/` itself holds only `__init__.py` and
`conftest.py`.

    uv run --extra test pytest tests/unit -q

161 tests, ~40 s. Never a bare `pytest`.

`conftest.py` chdirs each test module into a same-named subdirectory when one
exists (`tests/unit/test_rcsb/`, `tests/unit/test_pdbwrite/`). That matters
because `PDBParser.fetch()` downloads into the *current directory* and reuses a
file that is already there: the committed fixtures are what keep the suite
offline. A test naming an entry with no cached fixture will hit the network.

CI (`.github/workflows/tests.yaml`) runs pytest on 3.10/3.11/3.12 plus a
`sphinx -b doctest` job over `docs/source` — **doctest examples in the docs are
executed**, so an output change breaks CI even when the tests pass.

## The YAML format files are the spec

`pidibble/resources/pdb_format.yaml` (and `mmcif_format.yaml`) drive everything.
Adding or correcting a record type is normally a YAML edit, not a Python edit.

- A field spec is `[type, [start, end]]` with 1-based inclusive columns, plus an
  optional third element of **writer hints** (`{prec, just, sep}`) that the
  parser ignores. Prefer adding a hint over special-casing the writer.
- `PDBParser.record_formats` is the single active table that *both* parsing and
  writing read. Never let read and write diverge on columns — that invariant is
  what makes them exact inverses.
- Dialects: `_apply_dialect()` merges `charmm_formats` over the base table when
  `dialect='charmm'`. CHARMM widens `resName` to 6 columns and adds the segID
  column (73-76) while pinning x/y/z at 31-54; `ANISOU` is excluded from that
  dialect because its u-values would collide with the wide residue block.
- A local `pdb_format.yaml` in the CWD overrides the packaged one (used by
  `tests/unit/test_rcsb/test_pdb_format.yaml`).

## Round-trip contract

Not uniform across record families, by design:

- Single-line records (types 1/3) and `TER` round-trip **byte-exact**.
- Continuation/group records (types 2/4 — `COMPND`, `SEQRES`, `SITE`, …)
  round-trip **re-parse-identical, not byte-exact**: line breaks are discarded
  at parse time and the writer re-wraps. Don't chase byte-exactness there.
- Dispatch to multi-line emission on the presence of `continues` in the format,
  **not** on the type number (`REVDAT` is tagged type 3 but is multi-line).
- Hex serials flip decimal→hex permanently at the first serial >99999 and stay
  hex for the rest of the document, back-references included (plain hex, not
  hybrid-36); five serial columns puts the ceiling at `0xFFFFF` (1,048,575).
  Past it a serial is unrepresentable: VMD writes `*****`, the parser reads
  that marker as 0, and the writer emits it. **The marker test must stay ahead
  of the hex branch in `AtomSerialParser`** — hex trips at serial 100000, so a
  guard behind it can never fire on a real file — which is exactly the bug
  that made every structure over 1,048,575 atoms unreadable until 2026-09-20.
- From an mmCIF parse, `REMARK` and `COMPND`/`SOURCE` cannot be written — the
  parse does not retain what those records need. `REMARK`/`JRNL` are otherwise
  passed through verbatim from the source lines.

**When changing the writer, the correctness gate is the entry's own `.pdb` file**
from the RCSB, not the format spec: fetch both formats of one entry and diff the
emitted block against the real file. Several layout defects (separator style,
justification of fields whose declared type disagrees with their rendering) were
only ever caught that way.

## Citations

`PDBParser.citations()` / `.citation_ids()` **never open a socket**. `enrich=True`
is the only network path (RCSB Data API); a failed request is logged and
ignored, never raised — enrichment can only add. A test enforces this by making
`urlopen` raise. Keep it that way: a caller may be running inside an offline
build. PDB-sourced citations cannot recover author/title capitalization
(`A.B.MCDERMOTT`); titles are handed back as stored rather than guessed at.

## Releases

    scripts/release.sh <version>

It requires a clean tree on `main` with an `## [Unreleased]` section, rotates
the CHANGELOG, bumps `pyproject.toml`, commits, tags `v<version>`, and pushes.
**Pushing the tag is what publishes** (PyPI, GitHub Release, ReadTheDocs). Never
hand-roll any of those steps. Write changes under `## [Unreleased]` as you go,
Keep a Changelog style.

## Known rough edges

- **The floor is Python 3.10** (`requires-python`), and it is load-bearing:
  `pdbparse.py` uses PEP 604 unions (`str | Path`) in a signature without
  `from __future__ import annotations`, which evaluate at def time. Don't lower
  it — a too-low floor also breaks `uv run --extra test`, because uv resolves
  across the whole declared range and no `numpy>=1.24` supports 3.7-3.9. If a
  lower floor is ever wanted, the code has to change first.
- **`mmcif` is optional and must stay lazily imported.** Only
  `_require_mmcif()` in `pdbparse.py` may import it; a module-level import
  would make every install need a compiled package conda-forge lacks.
  `tests/unit/test_optional_mmcif.py` checks this in a fresh interpreter. The
  `test` extra pulls it in, so the suite still covers the mmCIF path.
- **Docs builds must install the checkout, never `pidibble` from PyPI.**
  `docs/source/conf.py` reads the version from installed metadata and autodoc
  imports the installed package. `.readthedocs.yaml` pip-installs `path: .` and
  CI's doctest job does `pip install -e .`; don't re-add `pidibble` to
  `docs/requirements.txt`.
- Docstrings are numpydoc; `docs/source/api/` is autodoc, so a signature change
  shows up in the published docs without any docs edit.
- `docs/testing-notes.md` explains why a broken serial round-trip makes the
  suite *appear* to hang (pytest's difflib rendering on 4500-element lists), not
  a parser loop.
