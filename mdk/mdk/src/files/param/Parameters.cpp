#include "files/param/Parameters.hpp"
using namespace mdk::param;
using namespace std;

Parameters::Parameters() {
    defAngleParams = Coeffs(7);
    for (auto const& var: pairTypes()) {
        angleParams[var] = defAngleParams;
        dihedralParams[var] = defAngleParams;
    }

    for (auto const& acid1: AminoAcid::aminoAcids()) {
        specificity[acid1] = SpecificityParams();
        radius[acid1] = 0;

        for (auto const& acid2: AminoAcid::aminoAcids()) {
            pairwiseMinDist[{acid1, acid2}] = 0;
        }
    }
}
