# Pidbibble
> Fully-functional PDB file parser

Pidibble is a Python package for parsing standard Protein Data Bank (PDB) files.  It conforms to the [https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html](most recent standard) (v.3.3 Atomic Coordinate Entry Format, ca. 2010).  Data from all allowed keywords is read in and parsed.

## Installation

```bash
git clone git@github.com:cameronabrams/pidibble.git
cd pidibble
pip install -e .
```

Once installed, the user has access to the `PDBParser` class in the `pidibble.rcsb` module.

## Release History

* 1.0-rc-1
    * Initial version

## Meta

Cameron F. Abrams – cfa22@drexel.edu

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

[https://github.com/AbramsGroup](https://github.com/AbramsGroup/)

## Contributing

1. Fork it (<https://github.com/AbramsGroup/HTPolyNet/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

