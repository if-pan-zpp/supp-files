#pragma once
#include "../../utils/AminoAcid.hpp"
#include "../cmap/ContactMap.hpp"
#include "../../model/Model.hpp"

namespace mdk::seq {
    /**
     * A sequence file, i.e. a set of chains with optional interposed contact
     * maps.
     */
    class Sequence {
    public:
        /// A single chain in the sequence file.
        struct Chain {
            /// Amino acids comprising the chain.
            std::vector<AminoAcid> codes;

            /**
             * Paths to contact maps associated with the chain. Technically
             * one could store \p cmap::ContactMap explicitly but we store paths
             * in order to not duplicate contact maps.
             */
            std::vector<std::string> contactMaps;
        };

        /// Optional screening distance, as per README.txt.
        std::optional<double> screeningDist;

        /// Chains in the sequence file.
        std::vector<Chain> chains;

        /// A map from the contact map path to a parsed contact map.
        std::unordered_map<std::string, cmap::ContactMap> contactMaps;

    public:
        /**
         * Create a \p Model from the sequence file. Positions are initialized
         * to zero, thus invoking a function like \p morphIntoSAW is necessary.
         * Contacts are assumed to be "NAT", i.e. without specifying which is
         * sidechain or backbone etc.
         * @return Resulting \p Model from the sequence file
         */
        Model asModel() const;
    };
}