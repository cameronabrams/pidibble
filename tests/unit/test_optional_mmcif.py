"""

.. module:: test_optional_mmcif
   :synopsis: pidibble works without the optional mmcif package, except for mmCIF input

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

``mmcif`` (wwPDB py-mmcif) is an optional extra. These tests hide it and check
that everything outside the mmCIF input path still works, and that asking for
mmCIF fails early and says how to fix it. Offline: they use the fixtures in
``test_rcsb/``.

"""
import subprocess
import sys

import pytest

from pidibble.pdbparse import PDBParser


@pytest.fixture
def no_mmcif(monkeypatch):
    """Make ``import mmcif...`` raise, as if the package were not installed."""
    for name in ('mmcif', 'mmcif.io', 'mmcif.io.IoAdapterCore'):
        monkeypatch.setitem(sys.modules, name, None)


def test_importing_pidibble_does_not_import_mmcif():
    # a fresh interpreter, so nothing else in the test session has imported it
    code = ('import sys, pidibble.pdbparse, pidibble.pdbwrite, pidibble.citation; '
            'print("mmcif" in sys.modules)')
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == 'False'


def test_pdb_parse_and_write_without_mmcif(no_mmcif):
    p = PDBParser(filepath='test_rcsb/4zmj.pdb').parse()
    assert len(p.parsed['ATOM']) == 4518
    assert p.write_PDB()[0].startswith('HEADER')
    assert p.citations()


def test_mmcif_input_without_mmcif_raises_before_fetch(no_mmcif, monkeypatch):
    def fetch_must_not_run(self):
        raise AssertionError('fetch() ran before the missing dependency was reported')
    monkeypatch.setattr(PDBParser, 'fetch', fetch_must_not_run)
    with pytest.raises(ImportError, match=r"pidibble\[mmcif\]"):
        PDBParser(source_db='rcsb', source_id='4zmj', input_format='mmCIF').parse()
