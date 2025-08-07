# Pidibble 
> a complete PDB-file parser

[![PyPI Downloads](https://static.pepy.tech/badge/pidibble)](https://pepy.tech/projects/pidibble)

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

# Example downloading from AlphaFold

```
>>> from pidibble.pdbparse import PDBParser
>>> p=PDBParser(alphafold='O46077').parse()
>>> p.parsed['TITLE'].title
'ALPHAFOLD MONOMER V2.0 PREDICTION FOR ODORANT RECEPTOR 2A (O46077)'
>>> print(p.parsed['ATOM'][0].pstr())
ATOM
              serial: 1
                name: N
              altLoc: 
             residue: resName: MET; chainID: A; seqNum: 1; iCode: 
                   x: -0.553
                   y: 26.513
                   z: 23.174
           occupancy: 1.0
          tempFactor: 39.74
             element: N
              charge: 


```

## Release History
* 1.4.1
   * Fixed parsing bug in `PDBRecordList`
* 1.4.0
   * Introduced `PDBRecordList` and `PDBRecordDict` classes
* 1.3.3
   * fixed bugs regarding assuming missing records actually present in mmcif parsing
* 1.3.2
   * implemented `filepath` parameter in `PDBParser()` to make
     reading local files more transparent
* 1.3.0
   * streamline class attribute usage; full API documentation
* 1.2.3:
   * bugfix: negative resids allowed
* 1.2.1:
   * bugfix: hex issues AGAIN
* 1.2.0
   * bugfix: hex issues again
* 1.1.9
   * gently handles any nan values or '*' fillers
* 1.1.8
   * Allow unstructured REMARK records (thanks a lot, packmol) to be parsed as REMARK.-1
* 1.1.6
   * bugfix: detection of hex when no abcdef present based on a trip
* 1.1.5
   * bugfix: allow for ONLY atom serial numbers to be hex-or-int
* 1.1.4
   * Added ability to read hexadecimal atom indices for files with > 99999 atoms
* 1.1.3
   * Added ability to group into models for multiple-model entries
* 1.1.2
   * Added ability to retrieve pdbs from AlphaFold database
* 1.1.1
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


