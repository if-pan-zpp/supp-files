#pragma once
#include <iostream>
#include "Parameters.hpp"

namespace mdk::param {
    /**
     * A parser of parameter files in the legacy format (as described in
     * legacy README.txt).
     */
    class LegacyParser {
    public:
        /**
         * Read parameter file from a stream.
         * @param is Input stream wherefrom to read the contact map
         * @return Parsed parameter file.
         */
        Parameters read(std::istream& is);

        /**
         * Write parameter file to a stream.
         * @param os Output stream whereto to write the contact map
         * @param data Parameter file to write
         * @return \p os
         */
        std::ostream& write(std::ostream& os, Parameters const& data);
    };
}