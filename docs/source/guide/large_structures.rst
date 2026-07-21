Large structures and hexadecimal serials
=========================================

The legacy PDB format allots five columns to the atom serial number, so it can
only count up to ``99999``.  Structures with more atoms than that continue the
count in **hexadecimal** (``100000`` is written ``186A0``), a widely used
convention that a naive integer parse would get wrong.

pidibble handles this automatically — there is nothing you need to configure.
Each :class:`~pidibble.pdbparse.PDBParser` owns an
:class:`~pidibble.hex.AtomSerialParser` that watches the atom serial column and
switches to hexadecimal parsing at the right moment, so ``serial`` (and the
serial references in ``CONECT``, ``TER`` and ``ANISOU``) stay correct past
``99999``:

.. code-block:: python

   p = PDBParser(source_db='rcsb', source_id='<large-entry>').parse()
   p.parsed['ATOM'][100000].serial      # a correct integer, not a mis-parsed hex string

How the switch is detected
--------------------------

The transition to hexadecimal is recognized two ways, so it works whether or not
the early serial numbers happen to contain the digits ``a``–``f``:

* **By content** — as soon as a serial contains a hexadecimal letter, subsequent
  serials are read as hexadecimal.
* **By magnitude** — if the running value exceeds ``99999``, the parser trips into
  hexadecimal mode even for an all-numeric field.

The switch is *latching*: once tripped it stays in hexadecimal mode for the rest
of the parse.  Because the state lives on the parser instance (not in a global),
parsing several structures in the same program — or in several threads — never
lets one file's numbering leak into another's.

.. note::

   Hexadecimal detection is applied **only** to atom-serial fields.  Ordinary
   integer fields (residue sequence numbers, counts, and the like) are always
   read as decimal, including legitimately negative values.
