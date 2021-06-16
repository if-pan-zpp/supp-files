#include "forces/Tether.hpp"
#include "data/Chains.hpp"
#include "simul/Simulation.hpp"
using namespace mdk;

Tether::Tether(bool fromNative) {
    this->fromNative = fromNative;
}

void Tether::bind(Simulation &simulation) {
    Force::bind(simulation);
    n = state->n;
    dist0 = Scalars::Constant(n, 3.8 * angstrom);
    isConnected = Bytes(n, false);

    auto model = simulation.data<Model>();
    for (auto const& chain: model.chains) {
        for (int i = chain.start; i + 1 < chain.end; ++i) {
            isConnected[i] = true;
        }

        if (fromNative) {
            auto const& residues = model.residues;
            for (auto i = chain.start; i + 1 < chain.end; ++i) {
                if (not residues[i].nat_r || not residues[i+1].nat_r) {
                    throw std::runtime_error("no native position"
                                             " - tether's length can't be calculated");
                }
                Vector r1 = *residues[i].nat_r, r2 = *residues[i+1].nat_r;
                auto r12 = r2 - r1;
                dist0[i] = r12.norm();
            }
        }
    }
}

void Tether::asyncPart(Dynamics &dyn) {
    #pragma omp for nowait
    for (int i = 0; i < n - 1; ++i) {
        if (not isConnected[i]) continue;
        
        auto r1 = state->r[i], r2 = state->r[i+1];
        auto r12 = r2 - r1;
        auto r12_norm = r12.norm();

        auto dx = r12_norm - dist0[i];
        auto r12_unit = r12 / r12_norm;
        harm.computeF(r12_unit, dx, dyn.V, dyn.F[i], dyn.F[i+1]);
    }
}
