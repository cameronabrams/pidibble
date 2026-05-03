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

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for the full release history.
