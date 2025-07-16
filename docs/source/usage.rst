Usage
=====

Usage Example
-------------

Let's parse the PDB entry '4ZMJ', which is a trimeric ectodomain construct of the HIV-1 envelope glycoprotein:

>>> from pidibble.pdbparse import PDBParser
>>> p=PDBParser(PDBcode='4zmj').parse()

The ``PDBParser()`` call creates a new ``PDBParser`` object, and the member function ``parse()`` executes (optionally) downloading the PDB file of the code entered with the ``PDBcode`` keyword argument to ``PDBParser()``, followed by parsing into a member dictionary ``parsed``.  (The file is downloaded from the RSCB only if it is not found in the current working directory.)

Alternatively, a ``PDBParser()`` invocation can fetch from the AlphaFold model database by providing the accession code with the ``alphafold`` keyword:

>>> p=PDBParser(alphafold='O46077').parse()

(Note that this fetches a model for the odor receptor OR2a from *D. melanogaster*.  For the rest of this example, we'll work with the HIV-1 Env trimer 4zmj above.)

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
               label: BIOMT, BIOMT, BIOMT
          coordinate: 1, 2, 3
           divnumber: 1, 1, 1
                 row: [m1: 1.0; m2: 0.0; m3: 0.0; t: 0.0], [m1: 0.0; m2: 1.0; m3: 0.0; t: 0.0], [m1: 0.0; m2: 0.0; m3: 1.0; t: 0.0]
              header: G, B, A, C, D
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
>>> M,T=get_symm_ops(b)
>>> print(str(M))
[[1. 0. 0.]
 [0. 1. 0.]
 [0. 0. 1.]]
>>> print(str(T))
[0. 0. 0.]
>>> b=p.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM2']
>>> M,T=get_symm_ops(b)
>>> print(str(M))
[[-0.5      -0.866025  0.      ]
 [ 0.866025 -0.5       0.      ]
 [ 0.        0.        1.      ]]
>>> print(str(T))
[107.18    185.64121   0.     ]
>>> b=p.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM3']
>>> M,T=get_symm_ops(b)
>>> print(str(M))
[[-0.5       0.866025  0.      ]
 [-0.866025 -0.5       0.      ]
 [ 0.        0.        1.      ]]
>>> print(str(T))
[-107.18     185.64121    0.     ]

You may recognize these rotation matrices as those that generate an object with C3v symmetry.  Each rotation is also accompanied by a translation, here in the ``Tlist`` object.

Because many entries in the RCSB do not have "legacy" PDB files and instead only have the (now standard) mmCIF/PDBx format files, ``pidibble`` can also generate parsed objects from these files.  This is activated by specifying a value ``mmCIF`` for the ``input_format`` keyword argument to the ``PDBParser`` generator:

>>> from pidibble.pdbparse import PDBParser
>>> p=PDBParser(PDBcode='4tvp',input_format='mmCIF').parse()
>>> b=p.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM2']
>>> print(b.pstr())
REMARK.350.BIOMOLECULE1.TRANSFORM2
         BIOMOLECULE: 1
           tmp_label: BIOMOLECULE1.TRANSFORM2
              header: A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T
           divnumber: 2
           TRANSFORM: 2
                row1: m1: -0.5; m2: -0.8660254038; m3: 0.0; t: -515.56
                row2: m1: 0.8660254038; m2: -0.5; m3: 0.0; t: 0.0
                row3: m1: 0.0; m2: 0.0; m3: 1.0; t: 0.0
                 row: [m1: -0.5; m2: -0.8660254038; m3: 0.0; t: -515.56], [m1: 0.8660254038; m2: -0.5; m3: 0.0; t: 0.0], [m1: 0.0; m2: 0.0; m3: 1.0; t: 0.0]
          coordinate: 1, 2, 3

We can compare this to the ``REMARK.350.BIOMOLECULE1.TRANSFORM2`` record from the analogous PDB file:

>>> p=PDBParser(PDBcode='4tvp',input_format='PDB').parse()
>>> b=p.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM2']
>>> print(b.pstr())
REMARK.350.BIOMOLECULE1.TRANSFORM2
               label: BIOMT, BIOMT, BIOMT
          coordinate: 1, 2, 3
           divnumber: 2, 2, 2
                 row: [m1: -0.5; m2: -0.866025; m3: 0.0; t: -515.56], [m1: 0.866025; m2: -0.5; m3: 0.0; t: 0.0], [m1: 0.0; m2: 0.0; m3: 1.0; t: 0.0]
              header: G, B, L, H, D, E, A, C, F, I, J, K, M, N, O, P, Q, R, S, T

Note that the important attributes of ``row`` and ``header`` are the same (in ``header``'s case, the lists are in different orders but they have the same elements).  Note the greater precision in the floating-point values for the record read in from the ``mmCIF`` file.

Currently, only ``ATOM``, ``HETATM``, ``SEQADV``, ``REMARK 350``, and ``REMARK 465`` records are translated from a ``mmCIF``-format file: 

>>> ', '.join(list(p.parsed.keys()))
'ATOM, HETATM, LINK, SSBOND, SEQADV, REMARK.350.BIOMOLECULE1.TRANSFORM1, REMARK.350.BIOMOLECULE1.TRANSFORM2, REMARK.350.BIOMOLECULE1.TRANSFORM3, REMARK.465'

These records are the bare minimum needed to generate (say) input coordinate and topology files for an MD simulation.  Future versions of ``pidibble`` will provide complete PDB-like parsings of ``mmCIF`` files.  This is probably not useful.

Importantly:  ``pidibble`` parses mmCIF input to generate a structure that is the equivalent of the PDB format; that is, it uses ``auth`` fields instead of ``label`` fields.  