#include "forces/dihedral/HeuresticDihedral.hpp"
#include "forces/dihedral/DihedralAngles.hpp"
using namespace mdk;

void HeuresticDihedral::bind(Simulation &simulation) {
    auto& model = simulation.data<Model>();
    auto& params = simulation.data<param::Parameters>();

    angleTypes = Eigen::Matrix<int8_t, Eigen::Dynamic, 1>(model.n);
    for (auto const& chain: model.chains) {
        auto start = chain.start + 2;
        auto end = chain.end - 1;

        for (int i = start; i < end; ++i) {
            auto i2 = i-1, i3 = i;
            AminoAcid acid2(int8_t(model.residues[i2].type));
            AminoAcid acid3(int8_t(model.residues[i3].type));

            PairType pt = pairType(acid2, acid3);
            angleTypes[i] = (int8_t)pt;
        }
    }

    for (auto const& [pt, coeffs]: params.dihedralParams) {
        for (int d = 0; d < 6; ++d) {
            coeff[(int8_t)pt][d] = coeffs[d];
        }
    }

    auto& unifiedDih = simulation.var<DihedralAngles>();
    unifiedDih.heurDih = this;
}
