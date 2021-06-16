#pragma once
#include <string>
#include <vector>
#include "AminoAcid.hpp"

namespace mdk {
    /**
     * A "pair type", as distinguished by a number of force fields, such
     * as \p HeuresticBondAngles or \p HeuresticDihedralAngles.
     */
    enum class PairType: int8_t {
        GG, GP, GX, PG, PP, PX, XG, XP, XX
    };

    constexpr inline int numOfPTs = 9;
    std::vector<PairType> pairTypes();
    PairType pairType(AminoAcid const& acid1, AminoAcid const& acid2);
}