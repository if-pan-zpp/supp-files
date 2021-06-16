#include "forces/dihedral/NativeDihedralBase.hpp"
#include "data/Chains.hpp"
using namespace mdk;

void NativeDihedralBase::bind(Simulation &simulation) {
    auto const& model = simulation.data<Model>();
    isNative = simulation.data<Chains>().nativeQuads;
    phi0 = Scalars(model.n);

    for (auto const& ch: model.chains) {
        for (auto const& spIdx: ch.structuredParts) {
            auto const& sp = model.structuredParts[spIdx];
            auto spStart = ch.start + sp.off + 2;
            auto spEnd = ch.start + sp.off + sp.len - 1;

            for (int i = spStart; i < spEnd; ++i) {
                phi0[i] = sp.dihedral[i - (ch.start + sp.off)];
            }
        }
    }
}