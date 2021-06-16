#include <filesystem>
#include <optional>
#include "../../utils/AminoAcid.hpp"
#include "Sequence.hpp"

namespace mdk::seq {
    /**
     * A parser of sequence files. Here we need an explicit path since the
     * sequence files reference contact maps by relative path.
     */
    class LegacyParser {
    public:
        /**
         * Read a sequence from a file at \p path.
         * @param path Path to a file with the sequence.
         * @return Parsed sequence.
         */
        Sequence read(std::filesystem::path const& path);
    };
}