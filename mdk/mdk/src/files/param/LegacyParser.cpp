#include "files/param/LegacyParser.hpp"
#include "utils/Text.hpp"
#include "utils/Units.hpp"
#include <sstream>
using namespace mdk;
using namespace mdk::param;
using namespace std;

static void fetchDefAngleParams(istream &is, Parameters &data) {
    skipLine(is);
    auto ss = lineStream(is);
    for (auto& coeff: data.defAngleParams)
        ss >> coeff;
}

static void fetchAngleParams(istream &is, Parameters &data) {
    skipLine(is);
    for (auto var: pairTypes()) {
        auto ss = lineStream(is);
        for (auto& coeff: data.angleParams[var])
            ss >> coeff;
    }
}

static void fetchDihedralParams(istream &is, Parameters &data) {
    skipLine(is);
    for (auto var: pairTypes()) {
        auto ss = lineStream(is);
        for (auto& coeff: data.dihedralParams[var])
            ss >> coeff;
    }
}

static void fetchSpecificity(istream &is, Parameters &data) {
    auto ss = lineStream(is);
    vector<AminoAcid> order;
    for (string name; ss >> name; ) {
        auto acid = (AminoAcid)(string)name;
        order.push_back(acid);
    }

    ss = lineStream(is);
    for (auto const& acid: order)
        ss >> (int8_t&)data.specificity[acid].polarization;

    ss = lineStream(is);
    for (auto const& acid: order)
        ss >> data.specificity[acid].maxSidechain;

    ss = lineStream(is);
    for (auto const& acid: order)
        ss >> data.specificity[acid].maxHydrophobicSS;

    ss = lineStream(is);
    for (auto const& acid: order)
        ss >> data.specificity[acid].maxPolarSS;
}

static void fetchRadii(istream &is, Parameters &data) {
    skipLine(is);

    auto ss = lineStream(is);
    for (auto const& acid: AminoAcid::aminoAcids())
        ss >> data.radius[acid];
}

static void fetchPairwiseData(istream &is, Parameters &data) {
    string header;
    getline(is, header);

    bool fetchMJ = (header == "amino acid pair distances and energies");
    if (fetchMJ) {
        data.mjMatrix = Parameters::PerPairData<double>();
    }

    int numPairs = AminoAcid::N * (AminoAcid::N + 1) / 2;
    for (int i = 0; i < numPairs; ++i) {
        auto ss = lineStream(is);

        string name1, name2;
        ss >> name1 >> name2;
        AminoAcid acid1(name1), acid2(name2);

        double minDist;
        ss >> minDist;
        minDist *= angstrom;
        data.pairwiseMinDist[{acid1, acid2}] = minDist;
        data.pairwiseMinDist[{acid2, acid1}] = minDist;

        if (fetchMJ) {
            double mjEnergy;
            ss >> mjEnergy;
            mjEnergy *= eps;

            data.mjMatrix.value()[{acid1, acid2}] = mjEnergy;
            data.mjMatrix.value()[{acid2, acid1}] = mjEnergy;
        }
    }
}

Parameters LegacyParser::read(istream &is) {
    Parameters data;

    fetchDefAngleParams(is, data);
    fetchAngleParams(is, data);
    fetchDihedralParams(is, data);
    fetchSpecificity(is, data);
    fetchRadii(is, data);
    fetchPairwiseData(is, data);

    return data;
}
