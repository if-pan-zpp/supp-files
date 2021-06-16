#pragma once
#include <memory>
#include "Field.hpp"
#include "files/pdb/Record.hpp"

namespace mdk::pdb {
    using namespace records;

    /**
     * A parser of a single record. Because of the "bidirectionality" of \p
     * Field structures and ability to swap out the pointers, we can elegantly
     * represent as a list of the fields in the line of a record, and the
     * parsing and writing of the record are essentially provided out-of-the-box.
     */
    class RecordParser {
    protected:
        /**
         * Underlying record, which serves as a base of variable (integers, strings
         * etc. of a record) to which field parsers write parsed values and from
         * which the field parsers read the values to write.
         */
        Record record;

    public:
        /**
         * A list of fields comprising the record.
         */
        std::vector<std::shared_ptr<Field>> fields;

        /**
         * Try to parse a string into a PDB record. If unsuccessful, a
         * \p std::monostate() is returned.
         * @param s String to parse into a PDB record.
         * @return Either a parsed PDB record or a \p std::monostate() if
         * the parsing was unsuccessful.
         */
        Record tryParse(std::string const& s) const;

        /**
         * Write out the record (in the text form) to an output stream.
         * @param os Output stream to write the PDB record to.
         * @param other PDB record to write to the stream.
         * @return \p os
         */
        std::ostream& write(std::ostream& os, Record const& other);
    };

    class RemarkParser: public RecordParser {
    public:
        RemarkParser();
    };

    class AtomParser: public RecordParser {
    public:
        AtomParser();
    };

    class HetatmParser: public RecordParser {
    public:
        HetatmParser();
    };

    class SSBondParser: public RecordParser {
    public:
        SSBondParser();
    };

    class Cryst1Parser: public RecordParser {
    public:
        Cryst1Parser();
    };

    class LinkParser: public RecordParser {
    public:
        LinkParser();
    };

    class ModelParser: public RecordParser {
    public:
        ModelParser();
    };

    class EndmdlParser: public RecordParser {
    public:
        EndmdlParser();
    };

    class EndParser: public RecordParser {
    public:
        EndParser();
    };

    class TerParser: public RecordParser {
    public:
        TerParser();
    };
}