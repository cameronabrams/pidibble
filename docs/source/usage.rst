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

>>> from pidibble.pdbparse import PDBParser
>>> p=PDBParser(PDBcode='4zmj').parse()
>>> print (p.parsed['HEADER'].classification)
VIRAL PROTEIN
>>> print (p.parsed['HEADER'].depDate)
04-MAY-15
>>> print (p.parsed['HEADER'].idCode)
4ZMJ
>>> keys=list(sorted(list(p.parsed.keys())))
>>> print(keys)
['ANISOU', 'ATOM', 'AUTHOR', 'CISPEP', 'COMPND', 'CONECT', 'CRYST1', 'DBREF', 'END', 'EXPDTA', 'FORMUL', 'HEADER', 'HELIX', 'HET', 'HETATM', 'HETNAM', 'JRNL.AUTH', 'JRNL.DOI', 'JRNL.PMID', 'JRNL.REF', 'JRNL.REFN', 'JRNL.TITL', 'KEYWDS', 'LINK', 'MASTER', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'REMARK.100', 'REMARK.2', 'REMARK.200', 'REMARK.280', 'REMARK.290', 'REMARK.290.CRYSTSYMMTRANS', 'REMARK.3', 'REMARK.300', 'REMARK.350', 'REMARK.350.BIOMOLECULE.1', 'REMARK.4', 'REMARK.465', 'REMARK.500', 'REVDAT', 'SCALE1', 'SCALE2', 'SCALE3', 'SEQADV', 'SEQRES', 'SHEET', 'SOURCE', 'SSBOND', 'TER', 'TITLE']
>>> header=p.parsed['HEADER']
>>> print(header.pstr())
HEADER
      classification: VIRAL PROTEIN
             depDate: 04-MAY-15
              idCode: 4ZMJ
>>> atoms=p.parsed['ATOM']
>>> len(atoms)
4518
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
>>> b=p.parsed['REMARK.350.BIOMOLECULE.1']
>>> print(b.pstr())
REMARK.350.BIOMOLECULE.1
             rowName: ['BIOMT2', 'BIOMT3', 'BIOMT1', 'BIOMT2', 'BIOMT3', 'BIOMT1', 'BIOMT2', 'BIOMT3']
             replNum: [1, 1, 2, 2, 2, 3, 3, 3]
                  m1: [0.0, 0.0, -0.5, 0.866025, 0.0, -0.5, -0.866025, 0.0]
                  m2: [1.0, 0.0, -0.866025, -0.5, 0.0, 0.866025, -0.5, 0.0]
                  m3: [0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
                   t: [0.0, 0.0, 107.18, 185.64121, 0.0, -107.18, 185.64121, 0.0]

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



.. To retrieve a list of random ingredients,
.. you can use the ``lumache.get_random_ingredients()`` function:

.. .. autofunction:: lumache.get_random_ingredients

.. The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
.. or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
.. will raise an exception.

.. .. autoexception:: lumache.InvalidKindError

.. For example:

.. >>> import lumache
.. >>> lumache.get_random_ingredients()
.. ['shells', 'gorgonzola', 'parsley']

