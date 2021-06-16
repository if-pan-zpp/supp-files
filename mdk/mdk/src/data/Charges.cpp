#include "data/Charges.hpp"
using namespace mdk;

Charges::Charges(const Model &model, param::Parameters const& params) {
    resize(model.n);

    for (int i = 0; i < model.n; ++i) {
        AminoAcid acid((int8_t)model.residues[i].type);
        auto pol = params.specificity.at(acid).polarization;

        if (pol == param::Polarization::POLAR_POS) {
            (*this)[i] = 1 * echarge;
        }
        else if (pol == param::Polarization::POLAR_NEG) {
            (*this)[i] = -1 * echarge;
        }
        else {
            (*this)[i] = 0 * echarge;
        }
    }
}