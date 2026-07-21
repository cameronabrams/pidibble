Nonconformant files and logging
================================

Real-world PDB files routinely deviate from the fixed-column standard — a
``charge`` column holding ``O-`` instead of a bare number, a record running past
the 80-column limit, a field that will not coerce to its declared type.  A single
file can contain tens of thousands of instances of the *same* deviation, so
pidibble does not log each one.  Instead it aggregates deviations by **type** and
prints a compact summary at the end of the parse.

The end-of-parse summary
------------------------

When a parse encounters nonconformances, pidibble emits one ``INFO`` line per
*type* — keyed by ``(record type, field, kind)`` — with a count and one concrete
exemplar, rather than one line per instance:

.. code-block:: text

   INFO  3 PDB-format nonconformance type(s) encountered while parsing 6cm3;
         set logging to DEBUG for per-instance detail:
   INFO    ATOM.charge: 24187 value(s) not coercible to Float — e.g. 'O-' at cols 79-80
   INFO    ...

Per-instance detail is still available — it is just demoted to ``DEBUG``.  Set
the logger to ``DEBUG`` to see every offending record:

.. code-block:: python

   import logging
   logging.basicConfig(level=logging.DEBUG)
   p = PDBParser(source_db='rcsb', source_id='6cm3').parse()

Inspecting nonconformances programmatically
-------------------------------------------

The registry is also available on the parser as
:attr:`~pidibble.pdbparse.PDBParser.nonconformances`, a
:class:`~pidibble.baseparsers.NonconformanceRegistry`.  It is truthy (and has a
non-zero ``len``) exactly when deviations were found, and
:meth:`~pidibble.baseparsers.NonconformanceRegistry.types` returns one tuple per
distinct type:

.. code-block:: python

   p = PDBParser(source_db='rcsb', source_id='6cm3').parse()

   if p.nonconformances:
       print(f'{len(p.nonconformances)} nonconformance type(s)')
       for record_key, field, kind, count, byte_range, example in p.nonconformances.types():
           print(f'  {record_key}.{field}: {count}× {kind} (e.g. {example!r})')

This makes it easy to fail a pipeline, warn a user, or record a data-quality
metric based on *what kinds* of nonconformance an entry has, without drowning in
per-atom noise.

.. note::

   The registry covers the **PDB** path.  On the mmCIF path, a coverage summary
   of *unmapped categories* is logged at parse time, but per-value coercion
   issues are not tracked in the registry.
