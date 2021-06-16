#include "utils/PairType.hpp"
#include "utils/AminoAcid.hpp"
using namespace mdk;
using namespace std;

vector<PairType> mdk::pairTypes() {
    return {
        PairType::GG, PairType::GP, PairType::GX,
        PairType::PG, PairType::PP, PairType::PX,
        PairType::XG, PairType::XP, PairType::XX
    };
}

static int8_t encode(AminoAcid const& acid) {
    if ((AminoAcidIdx)acid == AminoAcidIdx::GLY) return 0;
    else if ((AminoAcidIdx)acid == AminoAcidIdx::PRO) return 1;
    else return 2;
}

PairType mdk::pairType(const AminoAcid &acid1, const AminoAcid &acid2) {
    return (PairType)(3*encode(acid1)+encode(acid2));
}