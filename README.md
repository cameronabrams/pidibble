# Pidibble - a complete PDB parser

Pidibble is a Python package for parsing standard Protein Data Bank (PDB) files.  It conforms to the [most recent standard](https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html) (v.3.3 Atomic Coordinate Entry Format, ca. 2010).

Unlike parsers like that found in packages like [BioPython](https://biopython.org/wiki/PDBParser), `pidibble` provides meaningfully parsed objects from *all* standard PDB record types, not just ATOMs and CONECTs.

Once installed, the user has access to the `PDBParser` class in the `pidibble.rcsb` module.

## Release History

* 1.0
    * Initial version


