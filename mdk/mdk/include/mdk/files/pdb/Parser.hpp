#pragma once
#include <vector>
#include <memory>
#include <unordered_map>
#include "Data.hpp"

namespace mdk::pdb {
    /**
     * A forward decl to the \p RecordParser. It's in the private headers, but
     * for the \p std::shared_ptr it's enough to have a decl only.
     */
    class RecordParser;

    /// A parser of PDB files.
    class Parser {
    private:
        /// (Internal) set of record parsers.
        std::unordered_map<int, std::shared_ptr<RecordParser>> parsers;

    public:
        Parser();

        /**
         * Read PDB data from an input stream.
         * @param is Input stream wherefrom to read the data.
         * @return Read PDB data.
         */
        Data read(std::istream& is);

        /**
         * Write formatted PDB data to an output stream.
         * @param os Output stream whereto to write the PDB data.
         * @param data PDB data to write.
         * @return \p os.
         */
        std::ostream& write(std::ostream& os, Data const& data);
    };
}
