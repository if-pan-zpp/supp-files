PDB objects and parsing
=======================

Parsing and printing of PDB files is more complex than that of the other files.
There are following components:

- :cpp:class:`mdk::pdb::Parser` parses and prints :cpp:class:`mdk::pdb::Data`
  objects;
- :cpp:class:`mdk::pdb::Model` represents a full-atomic model as derived from
  the :cpp:class:`mdk::pdb::Data`, and allows for conversion from/to
  :cpp:class:`mdk::Model` and :cpp:class:`mdk::pdb::Data`
- :cpp:class:`mdk::pdb::Data` represents a parsed PDB file but in a raw format,
  i.e. not as a model;
- :cpp:type:`mdk::pdb::records::Record` is a single record in a PDB file, a
  variant type of many record types (see `records` for details).

In general, one takes an input stream, recovers :cpp:class:`mdk::pdb::Data` via
:cpp:class:`mdk::pdb::Parser`, converts it to a :cpp:class:`mdk::pdb::Model`,
then (after potentially modifying it) reduces it to a :cpp:class:`mdk::Model` to
pass to a :cpp:class:`mdk::Simulation`. If one wants to print
:cpp:class:`mdk::Model`, a :cpp:class:`mdk::pdb::Data` object can be generated
from one (via :cpp:class:`mdk::pdb::Model`), which can the written with the
parser.

We also provide a framework wherewith one can add custom records to parse,
see `parsing` for details.

.. toctree::
   :maxdepth: 2

   records
   raw
   model
   parsing
