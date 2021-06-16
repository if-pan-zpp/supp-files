#include "stats/Stats.hpp"
#include "simul/Simulation.hpp"
using namespace mdk;

void Stats::bind(Simulation &simulation) {
    auto& params = simulation.data<param::Parameters>();
    auto& model = simulation.data<Model>();

    stats = std::vector<Stat>(model.n);
    types = &simulation.data<Types>();
    polarization = std::vector<param::Polarization>(model.n);

    for (int i = 0; i < model.n; ++i) {
        auto& stat = stats[i];
        auto type = model.residues[i].type;
        auto& specificity = params.specificity.at(AminoAcid((int8_t)type));

        stat.backbone =  (ResTypeIdx)type != ResTypeIdx::PRO ? 2 : 1;
        stat.sidechain = specificity.maxSidechain;
        stat.hydrophobicSS = specificity.maxHydrophobicSS;
        stat.polarSS = specificity.maxPolarSS;
        polarization[i] = specificity.polarization;
    }
}
