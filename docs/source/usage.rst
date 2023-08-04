Usage
=====

.. _installation:

Installation
------------

To use Pidibble, first install it from PyPI:

.. code-block:: console

   (.venv) $ pip install pidibble

Usage Example
-------------

Let's parse the PDB entry '4ZMJ', which is a trimeric ectodomain construct of the HIV-1 envelope glycoprotein:
>>> from pidibble.pdbparse import PDBParser
>>> p=PDBParser(PDBcode='4zmj').parse()

We can easily ask what record types were parsed:
>>> list(sorted(list(p.parsed.keys())))
['ANISOU', 'ATOM', 'AUTHOR', 'CISPEP', 'COMPND', 'CONECT', 'CRYST1', 'DBREF', 'END', 'EXPDTA', 'FORMUL', 'HEADER', 'HELIX', 'HET', 'HETATM', 'HETNAM', 'JRNL.AUTH', 'JRNL.DOI', 'JRNL.PMID', 'JRNL.REF', 'JRNL.REFN', 'JRNL.TITL', 'KEYWDS', 'LINK', 'MASTER', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'REMARK.100', 'REMARK.2', 'REMARK.200', 'REMARK.280', 'REMARK.290', 'REMARK.290.CRYSTSYMMTRANS', 'REMARK.3', 'REMARK.300', 'REMARK.350', 'REMARK.350.BIOMOLECULE1.TRANSFORM1', 'REMARK.4', 'REMARK.465', 'REMARK.500', 'REVDAT', 'SCALE1', 'SCALE2', 'SCALE3', 'SEQADV', 'SEQRES', 'SHEET', 'SOURCE', 'SSBOND', 'TER', 'TITLE']

The `pstr()` method can be used to see the contents of any record:
>>> header=p.parsed['HEADER']
>>> print(header.pstr())
HEADER
      classification: VIRAL PROTEIN
             depDate: 04-MAY-15
              idCode: 4ZMJ

PDB records that are "single-line-multiple-occurrence", like ATOMs, HETATMs, SSBONDs, etc., are resolved as *lists* of pdbrecords:
>>> atoms=p.parsed['ATOM']
>>> len(atoms)
4518

Have a look at the first atom:
>>> print(atoms[0].pstr())
ATOM
              serial: 1
                name: N
              altLoc: 
             residue: resName: LEU; chainID: G; seqNum: 34; iCode: 
                   x: -0.092
                   y: 99.33
                   z: 57.967
           occupancy: 1.0
          tempFactor: 137.71
             element: N
              charge: 

Pidibble also parses any transformations needed to generate biological assemblies:
>>> b=p.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM1']
>>> print(b.pstr())
REMARK.350.BIOMOLECULE1.TRANSFORM1
             rowName: ['BIOMT2', 'BIOMT3', 'BIOMT1', 'BIOMT2', 'BIOMT3', 'BIOMT1', 'BIOMT2', 'BIOMT3']
             replNum: [1, 1, 2, 2, 2, 3, 3, 3]
                  m1: [0.0, 0.0, -0.5, 0.866025, 0.0, -0.5, -0.866025, 0.0]
                  m2: [1.0, 0.0, -0.866025, -0.5, 0.0, 0.866025, -0.5, 0.0]
                  m3: [0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
                   t: [0.0, 0.0, 107.18, 185.64121, 0.0, -107.18, 185.64121, 0.0]
              header: ['G', 'B', 'A', 'C', 'D']
              tokens:
AUTHOR DETERMINED BIOLOGICAL UNIT:  HEXAMERIC
SOFTWARE DETERMINED QUATERNARY STRUCTURE:  HEXAMERIC
            SOFTWARE USED:  PISA
TOTAL BURIED SURFACE AREA:  44090 ANGSTROM**2
SURFACE AREA OF THE COMPLEX:  82270 ANGSTROM**2
CHANGE IN SOLVENT FREE ENERGY:  81.0 KCAL/MOL

The `header` for any transform subrecord in a type-350 REMARK is the list of chains to which all transform(s) are
applied to generate this biological assembly.  If we send that record to the accessory method `get_symm_ops`, we can get `numpy.array()` versions of any matrices:
>>> from pidibble.pdbparse import get_symm_ops
>>> Mlist,Tlist=get_symm_ops(b)
>>> for m in Mlist:
...     print(str(m))
... 
[[1. 0. 0.]
 [0. 1. 0.]
 [0. 0. 1.]]
[[-0.5      -0.866025  0.      ]
 [ 0.866025 -0.5       0.      ]
 [ 0.        0.        1.      ]]
[[-0.5       0.866025  0.      ]
 [-0.866025 -0.5       0.      ]
 [ 0.        0.        1.      ]]

You may recognize these rotation matrices as those that generate an object C3v symmetry.  Each rotation is also accompanied by a translation, here in the `Tlist` object.