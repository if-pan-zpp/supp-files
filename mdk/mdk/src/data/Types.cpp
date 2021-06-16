#include "data/Types.hpp"
using namespace mdk;

Types::Types(const Model &model) {
    resize(model.n);
    for (int i = 0; i < model.n; ++i) {
        (*this)[i] = model.residues[i].type;
    }
}
