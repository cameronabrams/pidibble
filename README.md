# Pidibble - a complete PDB parser

Pidibble is a Python package for parsing standard Protein Data Bank (PDB) files.  It conforms to the [most recent standard](https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html) (v.3.3 Atomic Coordinate Entry Format, ca. 2011).

Unlike parsers like that found in packages like [BioPython](https://biopython.org/wiki/PDBParser), `pidibble` provides meaningfully parsed objects from *all* standard PDB record types, not just ATOMs and CONECTs.

Once installed, the user has access to the `PDBParser` class in the `pidibble.pdbparser` module.

# Example interactive usage

```
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

```
## Release History
* 1.1.0
   * version detection
* 1.0.9.1
   * added limited functionality to parse mmCIF files, in particular to generate any
     ATOM, HETATM, SSBOND, LINK, SEQADV, REMARK 350, and REMARK 465 records
* 1.0.8
    * bug fix: handle variations in how symmetry operation matrices are represented
* 1.0.7.7
    * cleaned up logging
* 1.0.7.6
    * bug fix: leading whitespace in resname field of Residue10 record sometimes ignored
* 1.0.7.5
    * support for four-letter residue names
* 1.0.7.4
    * added logging functionality
* 1.0.7.3
    * improved parsing of BIOMT transforms
* 1.0.7.2
    * added documentation stub at readthedocs
* 1.0.7.1
    * support for split BIOMT tables and REMARKS 280, 375, 650, and 700
* 1.0.7
    * pretty-print enabled
* 1.0
    * Initial version


