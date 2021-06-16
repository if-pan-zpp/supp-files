#pragma once
#include <vector>
#include <cstdint>
#include <string>

namespace mdk {
    /**
     * Underlying residue type index.
     */
    enum class ResTypeIdx: int8_t {
        ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY, HIS, ILE,
        LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL
    };

    /**
     * An object representing a residue type, along with a number of
     * facilities, conversions from/to other types etc.
     * At this moment, it is equivalent with the amino acids, but one may want
     * to add other pseudo-atoms, sites, water molecules or groups etc. so we
     * thought it wiser to separate the two.
     */
    class ResType {
    public:
        ResType() = default;

        constexpr explicit ResType(ResTypeIdx code): code {code} {};
        explicit operator ResTypeIdx const&() const;

        constexpr ResType(int8_t x): ResType((ResTypeIdx)x) {};
        operator int8_t() const;

        explicit ResType(std::string const& name);
        explicit operator std::string() const;

        double mass() const;

    private:
        ResTypeIdx code;
        friend struct std::hash<ResType>;
    };
}

namespace std {
    /**
     * A \p std::hash instantiation for \p ResType, allows us to use
     * \p ResType as keys in STL maps.
     */
    template<>
    struct hash<mdk::ResType> {
        size_t operator()(mdk::ResType const &resType) const {
            return std::hash<mdk::ResTypeIdx>()(resType.code);
        }
    };
}