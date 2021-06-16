#include "data/Masses.hpp"
using namespace mdk;

Masses::Masses(const Model &model) {
    resize(model.n);
    for (int i = 0; i < model.n; ++i) {
        (*this)[i] = model.residues[i].mass;
    }
}
