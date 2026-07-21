Advanced usage and customization
================================

pidibble is driven by a declarative description of the PDB and mmCIF formats, and
several constructor arguments let you extend or override it.

Custom field mappers
--------------------

Every field in the format specification names a *type*, and each type maps to a
callable that turns the raw column text into a Python value.  The built-in types
are ``String``, ``Integer``, ``Float`` and ``HxInteger`` (the hexadecimal-aware
atom-serial type), plus a family of list parsers.

You can add or override these with the ``mappers`` argument — a dict from type
name to a one-argument callable.  This is the hook for a custom type used by a
custom format, or for changing how an existing type is coerced:

.. code-block:: python

   def tenths(text):
       return round(float(text), 1)

   p = PDBParser(source_db='rcsb', source_id='4zmj',
                 mappers={'Float': tenths}).parse()

Your mappers are merged over the defaults, so you only supply what you want to
change.

Comment characters
------------------

By default, lines beginning with ``#`` are treated as comments and skipped.  Pass
``comment_chars`` to change that set — for example, to also ignore ``!`` lines:

.. code-block:: python

   p = PDBParser(filepath='annotated.pdb', comment_chars=['#', '!']).parse()

Overriding the format specification
-----------------------------------

The PDB and mmCIF record definitions are shipped as YAML resources
(``pdb_format.yaml`` and ``mmcif_format.yaml``).  To parse a nonstandard dialect
or add a record type, point the parser at your own files:

.. code-block:: python

   p = PDBParser(filepath='custom.pdb',
                 pdb_format_file='my_pdb_format.yaml').parse()

   c = PDBParser(filepath='custom.cif', input_format='mmCIF',
                 mmcif_format_file='my_mmcif_format.yaml').parse()

The format-specification schema
-------------------------------

Each record type in ``pdb_format.yaml`` is described by a small set of keys:

``type``
   The record's cardinality/shape: one-time-one-line, one-time-multiple-lines,
   multiple-times-one-line, multiple-times-multiple-lines, grouping, or other.

``fields``
   A mapping of field name to a ``[type, [start, end]]`` pair giving the field's
   data type and its (1-based, inclusive) column range.

``continues``
   The fields that a continuation record appends to or extends.

``token_formats``
   How a field's text is broken into named tokens (with optional ``determinants``
   that group token/value pairs) — this is what produces the token groups of
   ``COMPND`` and ``SOURCE``.

``concatenate``
   New fields built by list-concatenating other fields.

``allowed``
   Per-field lists of allowed values, for validation.

``determinants``
   The fields whose values decide whether a line starts a new record instance or
   continues the current one.

``subrecords``
   How to parse variant sub-formats, selected by the value of a ``branchon``
   field.

``tables``
   How embedded tables (such as the ``BIOMT`` matrices) are parsed out.

Reusable sub-formats (for example the ``Residue`` variants and the ``Biomt`` row)
live under ``custom_formats`` and are referenced by name from field definitions.
