# Pidibble 
> a complete PDB-file parser

[![PyPI Downloads](https://static.pepy.tech/badge/pidibble)](https://pepy.tech/projects/pidibble)

Pidibble is a Python package for parsing standard Protein Data Bank (PDB) files.  It conforms to the [most recent standard](https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html) (v.3.3 Atomic Coordinate Entry Format, ca. 2011).

Unlike parsers like that found in packages like [BioPython](https://biopython.org/wiki/PDBParser), `pidibble` provides meaningfully parsed objects from *all* standard PDB record types, not just ATOMs and CONECTs.

Once installed, the user has access to the `PDBParser` class in the `pidibble.pdbparse` module.

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
>>> keys=sorted(p.parsed.keys())
>>> len(keys)
60
>>> keys[:8]
['ANISOU', 'ATOM', 'AUTHOR', 'CISPEP', 'COMPND', 'CONECT', 'CRYST1', 'DBREF']
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
>>> p=PDBParser(source_db='alphafold', source_id='O46077').parse()
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
