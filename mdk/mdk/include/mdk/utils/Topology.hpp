#pragma once
#include "../data/Primitives.hpp"

namespace mdk {
    /**
     * An object representing the "topology" of the simulation box.
     * In particular, we want to compute PBC-adjusted distances between
     * residues.
     */
    class Topology {
    public:
        Vector cell, cellInv;
        bool use[3] = { false, false, false };

        void setCell(VRef cell);

        inline void fix(Vector& v) const {
            for (int i = 0; i < 3; ++i) {
                if (use[i])
                    v[i] -= round(v[i] * cellInv[i]) * cell[i];
            }
        }

        inline Vector operator()(Vector v) const {
            fix(v);
            return v;
        }
    };
}
