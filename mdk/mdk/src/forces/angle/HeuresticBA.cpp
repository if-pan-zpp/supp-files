#include "forces/angle/HeuresticBA.hpp"
#include "simul/Simulation.hpp"
#include "forces/angle/BondAngles.hpp"
using namespace mdk;

void HeuresticBA::bind(Simulation &simulation) {
    auto& model = simulation.data<Model>();
    auto& params = simulation.data<param::Parameters>();

    angleTypes = Bytes(model.n, 0);
    for (auto const& chain: model.chains) {
        for (int i = chain.start+1; i+1 < chain.end; ++i) {
            auto i2 = i, i3 = i+1;
            AminoAcid acid2(int8_t(model.residues[i2].type));
            AminoAcid acid3(int8_t(model.residues[i3].type));

            PairType pt = pairType(acid2, acid3);
            angleTypes[i] = (int8_t)pt;
        }
    }

    for (auto const& [pt, coeffs]: params.angleParams) {
        for (int d = 0; d < 7; ++d) {
            coeff[(int8_t)pt][d] = coeffs[d];
        }
    }

    auto& unifiedBA = simulation.var<BondAngles>();
    unifiedBA.heurBA = this;
}
