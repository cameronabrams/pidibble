# pidibble roadmap

Where an idea about this code waits before anyone acts on it. This is the
**only** place such items live: the fleet-level roadmap
(`~/.config/fleet/ROADMAP.md`) indexes this file, it never copies from it.

Nothing here is a commitment or a schedule. An item earns a place by being a
real, stated gap — not by being imaginable. Items move out by being done (note
the release) or by being dropped (say why).

Created 2026-10-04, at release 1.12.2.

## Open

### mmCIF coverage

- **`hydrog` and `saltbr` are dropped from `struct_conn`.** LINK folds in
  `covale`, `disulf` and `metalc` (2026-07-14); these two still fall on the
  floor. Unclear whether PDB LINK is the right target for either — decide that
  before mapping. Source: `docs/mmcif_coverage.md` item 6.
- **Evaluate `mmcif.api.DictionaryApi`** for type-correct coercion and for
  generating the format map rather than hand-maintaining it. Source:
  `docs/mmcif_coverage.md` item 11, already filed there as "later".

### Writer

- **`REMARK` and `COMPND`/`SOURCE` cannot be written from an mmCIF parse** —
  the parse does not retain what those records need. Closing this means
  deciding what the mmCIF path should retain, not patching the writer. Source:
  `CLAUDE.md`, round-trip contract.

## Done

*(nothing yet — items land here with the release that carried them)*

## Dropped

*(nothing yet — items land here with the reason)*
