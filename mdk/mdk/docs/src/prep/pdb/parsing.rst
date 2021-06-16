PDB parsing framework
=====================

The design of the PDB parsing framework is as follows:

- The field parsers provide the means of parsing from and writing to a PDB line
  of a single field (say, an integer or a real number with given width and
  precision etc.);
- A single PDB record parser is an ensemble of field parsers;
- A PDB parser has a collection of PDB record parsers which it tries to apply
  to a line in succession until it succeeds.

Adding a new record type
------------------------

In order to add a new record type, one must:

- (Potentially) add new field types;
- Add a record data structure to the list of records. (See `records` for
  details);
- Add a new derived class of :cpp:class:`mdk::pdb::RecordParser` for a given
  record type;
- Implement the constructor of the record parser, which in particular should
  populate the list of fields;
- Add the parser to the :cpp:class:`mdk::pdb::Parser` (in the constructor).

.. toctree::
   :maxdepth: 2
   :caption: Modules

   fields
   recordparsers