#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace mdk {
    /**
     * An underlying amino acid index.
     */
    enum class AminoAcidIdx: int8_t {
        ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY, HIS, ILE,
        LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL
    };

    /**
     * Info pertaining to a single atom type in an amino acid (for example
     * SG in cysteine).
     */
    struct AAAtomInfo {
        double radius;
        bool inBackbone = false;
    };

    /**
     * Info pertaining to an amino acid type.
     */
    struct AminoAcidInfo {
        double mass;

        /**
         * A set of names of the heavy atoms for the amino acid.
         */
        std::unordered_set<std::string> heavyAtoms;

        /**
         * A map from the heavy atom names to the atom info structures.
         */
        std::unordered_map<std::string, AAAtomInfo> atomInfo;
    };

    /**
     * An object representing an amino acid type, along with a number of
     * facilities, conversion from/to other types, general static data etc.
     */
    class AminoAcid {
    public:
        AminoAcid() = default;
        static std::vector<AminoAcid> aminoAcids();

        static std::vector<AminoAcidIdx> types();
        constexpr explicit AminoAcid(AminoAcidIdx type): type{type} {};
        explicit operator AminoAcidIdx const&() const;

        constexpr AminoAcid(int8_t x): AminoAcid((AminoAcidIdx)x) {};
        operator int8_t() const;

        explicit AminoAcid(char code);
        static std::string codes();
        explicit operator char const&() const;

        /**
         * Check whether a string corresponds to an amino acid proper name.
         * @param name Text to check.
         * @return true if \p name corresponds to an amino acid, false
         * otherwise.
         */
        static bool isProper(std::string const& name);
        explicit AminoAcid(std::string const& name);
        static std::vector<std::string> names();
        explicit operator std::string() const;

        bool operator==(AminoAcid const& other) const;
        bool operator!=(AminoAcid const& other) const;
        bool operator<(AminoAcid const& other) const;
        bool operator<=(AminoAcid const& other) const;
        bool operator>(AminoAcid const& other) const;
        bool operator>=(AminoAcid const& other) const;

        AminoAcidInfo const& info() const;

        static constexpr const int N = 20;

    private:
        AminoAcidIdx type;
        friend struct std::hash<AminoAcid>;
    };
}

namespace std {
    /**
     * A \p std::hash instantiation for \p AminoAcid, allows us to use
     * \p AminoAcid as keys in STL maps.
     */
    template<>
    struct hash<mdk::AminoAcid> {
        size_t operator()(mdk::AminoAcid const &aminoAcid) const {
            return std::hash<mdk::AminoAcidIdx>()(aminoAcid.type);
        }
    };
}
