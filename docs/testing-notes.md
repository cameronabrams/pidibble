# Testing notes

**Status:** note as of 2026-07-31 (pidibble 1.8.0). Developer reference, not
end-user documentation. Everything below has since been acted on; see
"Resolution" under each item.

## A failing large-list assertion makes the suite appear to hang

### Symptom

If a regression breaks serial-number round-tripping, running the unit suite
appears to hang indefinitely rather than report a failure:

```
tests/unit/test_pdbwrite.py::test_hex_serial_encoder_mirrors_parser PASSED [ 61%]
tests/unit/test_pdbwrite.py::test_big_serial_document_roundtrip
     <no further output; one core pinned at 100%, runs past 15 minutes>
```

The process is CPU-bound, not blocked on I/O or the network, which makes it look
like an infinite loop in the parser.

### Cause

It is not a pidibble bug. The cost is entirely inside pytest's assertion-diff
rendering.

`test_big_serial_document_roundtrip` ends with whole-list comparisons:

```python
assert [r.serial for r in parser.parsed['ATOM']] == [r.serial for r in q.parsed['ATOM']]
```

For 4ZMJ that is **4518 elements per side**. While the lists are equal the
assertion is cheap. Once a regression makes them differ, pytest builds a
human-readable diff through `difflib`, and `difflib._fancy_replace` recurses via
`_fancy_helper` computing pairwise similarity ratios across the mismatched
block. On lists this large, made of similar-looking values (`100001`, `100002`,
…), that descent is effectively unbounded.

Observed stack at the point of the "hang":

```
difflib.py:642  in quick_ratio
difflib.py:938  in _fancy_replace
difflib.py:997  in _fancy_helper
difflib.py:985  in _fancy_replace          <- recursing
...
_pytest/assertion/_compare_sequence.py:33  in _compare_eq_iterable
tests/unit/test_pdbwrite.py:262            in test_big_serial_document_roundtrip
```

### Confirmation

Disabling assertion introspection removes the cost completely, while the test
still fails as it should:

```console
$ pytest tests/unit/test_pdbwrite.py::test_big_serial_document_roundtrip -q
     <no completion; killed after 120s>

$ pytest tests/unit/test_pdbwrite.py::test_big_serial_document_roundtrip -q --assert=plain
1 failed in 0.72s
```

That is the diagnostic to reach for first: if `--assert=plain` returns
instantly, the time is going into diff rendering, not into pidibble.

Note that `-x` also hides the problem, because earlier tests in the file fail
first and stop the run before this test is reached. A regression can therefore
look like an ordinary failure under `-x` and like a hang without it.

### Resolution (applied)

`test_pdbwrite.py` now routes all three whole-list serial comparisons
(`ATOM`, `TER`, `CONECT`) through a `_assert_serials_equal(expected, got, label)`
helper that checks the length, then collects mismatches itself and reports only
a bounded slice — so `difflib` is never handed a large similar-element diff:

```python
assert len(got) == len(expected), f'{label}: got {len(got)} records, expected {len(expected)}'
mismatch = [(i, e, g) for i, (e, g) in enumerate(zip(expected, got)) if e != g]
assert not mismatch, f'{label}: {len(mismatch)} mismatches; first 10: {mismatch[:10]}'
```

The test still fails for exactly the same regressions. Re-checked against the
column mutation below: it now fails in **0.55s** instead of hanging, reporting

```
AssertionError: ATOM: 4518 mismatches; first 10: [(0, 100001, 34465), (1, 100002, 34466), ...]
```

which points straight at the first bad serial.

### How this was found

While using the suite as a measurement target, one-line mutations were applied
to source files to generate realistic red runs. The mutation that exposed this
was in the fixed-column field slice in `baseparsers.StringParser.parse`:

```python
# original (columns are 1-based in the PDB spec)
fieldstring = record[byte_range[0] - 1:byte_range[1]]
# mutated: drops the 1-based correction, shifting every field by one column
fieldstring = record[byte_range[0]:byte_range[1]]
```

Any regression that garbles parsed serials will reproduce it; that particular
mutation is just a convenient trigger.

Two unrelated observations from the same exercise, with what came of them:

- Two plausible mutations **survived** the suite — `safe_float`'s `'nan'` guard
  in `baseparsers.py` could be deleted, and the `> 99999` hex-trip threshold in
  `hex.AtomSerialParser.__call__` could be changed to `>= 99999`, with all 103
  tests still passing. Both were boundary/sentinel paths with no test.

  **Resolution (applied).** Both are now covered, and writing the tests turned
  up a real defect behind the first one. `safe_float` compared the *raw* field
  string to `'nan'`, but fixed-column float fields are right-justified, so a
  NaN written by an upstream tool arrives as `'     nan'` — the guard almost
  never fired, and `float('nan')` then let a NaN propagate silently into
  coordinates and B-factors. The comparison now strips and case-folds first
  (`x.strip().lower() == 'nan'`). Covered by
  `tests/unit/test_baseparsers.py` (direct, plus through `StringParser.parse`
  on a record with `nan` in a coordinate column) and by
  `test_hex_trip_threshold_is_strictly_above_99999` in `tests/unit/test_hex.py`,
  which pins that parsing `99999` must *not* trip hex mode — under `>=`, a
  following `CONECT` back-reference of `'10'` would silently read as 16. Each
  mutation was re-applied to confirm the new tests fail on it.

- The network-dependence claim recorded here originally was **backwards**, and
  is corrected as follows. `test_rcsb.py` (72 tests) and `test_hex.py` run
  fully offline: every structure they name — `1ca2`, `4tvp`, `4zmj`,
  `4zmj-newresnames`, `6m0j`, `8fae`, `test`, `my_system` — is a committed
  fixture under `tests/unit/test_rcsb/` or `tests/unit/test_hex/`, and
  `PDBParser.fetch()` only downloads when the file is absent. The one module
  that *did* hit RCSB was `test_pdbwrite.py`, whose `4zmj.pdb` was gitignored,
  so the byte-exactness round-trip ran against the live entry.

  **Resolution (applied).** That fixture is now committed under
  `tests/unit/test_pdbwrite/` and removed from the directory's `.gitignore`.
  The whole unit suite therefore runs offline on a clean checkout, and the
  byte-exact assertions are pinned to a known release rather than to whatever
  RCSB currently serves — a re-release of 4ZMJ can no longer break them for
  reasons unrelated to the writer. To pick up a revised entry deliberately,
  delete the file and re-run the suite; `fetch()` will download it again.
