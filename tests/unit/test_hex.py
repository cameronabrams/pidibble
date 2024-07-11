import numpy as np
from pidibble.pdbparse import PDBParser, PDBRecord, get_symm_ops
from pidibble.hex import int_or_hex
import unittest
import logging

def test_hex():
    a=int_or_hex('186a0')
    assert a==100000

def test_hex_in_pdb():
    p=PDBParser(PDBcode='my_system').parse()
    atoms=p.parsed['ATOM']
    assert len(atoms)==6
    a1=atoms[0]
    assert a1.serial==10000
    a1=atoms[1]
    assert a1.serial==99999
    a1=atoms[2]
    assert a1.serial==440691
