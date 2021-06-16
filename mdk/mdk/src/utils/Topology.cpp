#include "utils/Topology.hpp"
using namespace mdk;

void Topology::setCell(VRef _cell) {
    this->cell = _cell;
    for (int i = 0; i < 3; ++i) {
        if (cell[i] != 0.0) {
            cellInv[i] = 1.0 / cell[i];
        }
    }
}
