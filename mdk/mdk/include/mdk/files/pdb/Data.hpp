#pragma once
#include <vector>
#include "Record.hpp"
#include "Model.hpp"

namespace mdk::pdb {
    /**
     * PDB file in the "raw format", i.e. as a list of PDB records. Certain
     * utility functions, like parsing from/to \p pdb::Model and retrieving only
     * atom positions, are provided.
     */
    class Data {
    public:
        /// List of PDB records in the PDB file.
        std::vector<records::Record> records;

    public:
        Data() = default;

        /**
         * Convert a PDB model into a set of records. In particular, extracts
         * atom positions into ATOM records, contacts into SSBOND and LINK
         * records, cell shape into CRYST1.
         * @param model PDB model to convert.
         */
        explicit Data(pdb::Model const& model);

        /**
         * Extract only ATOM records from a PDB model. Useful for outputting
         * a single file with lots of conformations of the same amino acid
         * chain, for example when we output PDB model with consecutive
         * positions in the simulation.
         * @param model PDB model to convert.
         * @return PDB data object containing only ATOM records.
         */
        static Data onlyAtoms(pdb::Model const& model);

        /**
         * Convert the data into a PDB model (ATOM records into atoms,
         * SSBOND and LINK records into contacts, CRYST1 into cell shape).
         * @return PDB model generated from the data.
         */
        pdb::Model asModel() const;

        /**
         * Add a record to the PDB data object.
         * @param data PDB data object.
         * @param record PDB record to add.
         * @return \p data but with added \p record
         */
        friend Data& operator<<(Data& data, records::Record const& record);

        /**
         * Append a list of records from another PDB data object.
         * @param data PDB data object
         * @param other Another PDB data object from where to get the records
         * @return \p data but with added records from \p other
         */
        friend Data& operator<<(Data& data, Data const& other);
    };
}
