import numpy as np
from pidibble.pdbparse import PDBParser, PDBRecord, get_symm_ops
from pidibble.hex import str2atomSerial, AtomSerialParser
import unittest
import logging

def test_hex():
    a=str2atomSerial('186a0')
    assert a==100000

def test_hex_trip_threshold_is_strictly_above_99999():
    """99999 is the largest serial the 5-column field holds in decimal, so
    parsing it must *not* trip the parser into hex mode; only a value past the
    ceiling may. Getting the boundary wrong silently reinterprets every
    subsequent field -- a CONECT back-reference of '10' would become 16."""
    p=AtomSerialParser()
    assert p('99999')==99999
    assert not p._hex_tripped
    assert p('10')==10          # still decimal after the boundary value
    # one past the ceiling trips the parser, permanently
    assert p('100000')==100000
    assert p._hex_tripped
    assert p('10')==16          # every later field is now read as hex

def test_hex_in_pdb():
    p=PDBParser(PDBcode='my_system').parse()
    atoms=p.parsed['ATOM']
    assert len(atoms)==7
    a1=atoms[0]
    assert a1.serial==10000
    a1=atoms[1]
    assert a1.serial==99999
    a1=atoms[2]
    assert a1.serial==440691
    a1=atoms[3]
    assert a1.serial==65536
