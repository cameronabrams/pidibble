Usage
=====

.. _installation:

Installation
------------

To use Pidibble, install it from PyPI:

.. code-block:: console

   (.venv) $ pip install pidibble

Pidibble is also available via ``conda``: 

.. code-block:: console

   (conda-env) $ conda install -c conda-forge pidibble

Usage Example
-------------

Let's parse the PDB entry '4ZMJ', which is a trimeric ectodomain construct of the HIV-1 envelope glycoprotein:

>>> from pidibble.pdbparse import PDBParser
>>> p=PDBParser(PDBcode='4zmj').parse()

The ``PDBParser()`` call creates a new ``PDBParser`` object, and the member function ``parse()`` executes (optionally) downloading the PDB file of the code entered with the ``PDBcode`` keyword argument to ``PDBParser()``, followed by parsing into a member dictionary ``parsed``.

>>> type(p.parsed)
<class 'dict'>

We can easily ask what record types were parsed:

>>> list(sorted(list(p.parsed.keys())))
['ANISOU', 'ATOM', 'AUTHOR', 'CISPEP', 'COMPND', 'CONECT', 'CRYST1', 'DBREF', 'END', 'EXPDTA', 'FORMUL', 'HEADER', 'HELIX', 'HET', 'HETATM', 'HETNAM', 'JRNL.AUTH', 'JRNL.DOI', 'JRNL.PMID', 'JRNL.REF', 'JRNL.REFN', 'JRNL.TITL', 'KEYWDS', 'LINK', 'MASTER', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'REMARK.100', 'REMARK.2', 'REMARK.200', 'REMARK.280', 'REMARK.290', 'REMARK.290.CRYSTSYMMTRANS', 'REMARK.3', 'REMARK.300', 'REMARK.350', 'REMARK.350.BIOMOLECULE1.TRANSFORM1', 'REMARK.4', 'REMARK.465', 'REMARK.500', 'REVDAT', 'SCALE1', 'SCALE2', 'SCALE3', 'SEQADV', 'SEQRES', 'SHEET', 'SOURCE', 'SSBOND', 'TER', 'TITLE']

Every value in ``p.parsed[]`` is either a single instance of the class ``PDBRecord`` or a *list* of ``PDBRecords``.  Let's see which ones are lists:

>>> [x for x,v in p.parsed.items() if type(v)==list]
['REVDAT', 'DBREF', 'SEQADV', 'SEQRES', 'HET', 'HETNAM', 'FORMUL', 'HELIX', 'SHEET', 'SSBOND', 'LINK', 'CISPEP', 'ATOM', 'ANISOU', 'TER', 'HETATM', 'CONECT']

These are the so-called *multiple-entry* records; conceptually, they signify objects that appear more than once in a structure or it metadata.  Other keys each have only a single ``PDBRecord`` instance:

>>> [x for x,v in p.parsed.items() if type(v)!=list] 
['HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS', 'EXPDTA', 'AUTHOR', 'JRNL.AUTH', 'JRNL.TITL', 'JRNL.REF', 'JRNL.REFN', 'JRNL.PMID', 'JRNL.DOI', 'REMARK.2', 'REMARK.3', 'REMARK.4', 'REMARK.100', 'REMARK.200', 'REMARK.280', 'REMARK.290', 'REMARK.300', 'REMARK.350', 'REMARK.465', 'REMARK.500', 'CRYST1', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3', 'MASTER', 'END', 'REMARK.290.CRYSTSYMMTRANS', 'REMARK.350.BIOMOLECULE1.TRANSFORM1']
>>> type(p.parsed['HEADER'])
<class 'pidibble.pdbrecord.PDBRecord'>
>>> 

To get a feeling for what is in each record, use the ``pstr()`` method on any ``PDBRecord`` instance: 

>>> header=p.parsed['HEADER']
>>> print(header.pstr())
HEADER
      classification: VIRAL PROTEIN
             depDate: 04-MAY-15
              idCode: 4ZMJ

The format of this output tells you the instance attributes and their values:

>>> header.classification
'VIRAL PROTEIN'
>>> header.depDate
'04-MAY-15'
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

The ``header`` instance attribute for any transform subrecord in a type-350 REMARK is the list of chains to which all transform(s) are
applied to generate this biological assembly.  If we send that record to the accessory method ``get_symm_ops()``, we can get ``numpy.array()`` versions of any matrices:

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

You may recognize these rotation matrices as those that generate an object C3v symmetry.  Each rotation is also accompanied by a translation, here in the ``Tlist`` object.