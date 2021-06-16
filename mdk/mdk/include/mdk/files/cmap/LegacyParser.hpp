#pragma once
#include <iostream>
#include "ContactMap.hpp"

namespace mdk::cmap {
    /**
     * A parser of contact map files in the legacy format (as described in
     * legacy README.txt).
     */
    class LegacyParser {
    public:
        /**
         * Read contact map from a stream.
         * @param is Input stream wherefrom to read the contact map
         * @return Parsed contact map.
         */
        ContactMap read(std::istream& is);

        /**
         * Write contact map to a stream.
         * @param os Output stream whereto to write the contact map
         * @param cmap Contact map to write
         * @return \p os
         */
        std::ostream& write(std::ostream& os, ContactMap const& cmap);
    };
}
