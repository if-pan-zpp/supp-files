#include "files/seq/LegacyParser.hpp"
#include "files/cmap/LegacyParser.hpp"
#include <fstream>
using namespace mdk::seq;
using namespace std;

Sequence LegacyParser::read(const std::filesystem::path &path) {
    auto file = ifstream(path);
    Sequence seq;
    cmap::LegacyParser cmapParser;

    string line;
    getline(file, line);

    int numOfChains;
    if (line.substr(0, 7) == "screend") {
        seq.screeningDist = stod(line.substr(8));
        file >> numOfChains;
    }
    else {
        numOfChains = stoi(line);
    }

    seq.chains = vector<Sequence::Chain>(numOfChains);
    for (auto& chain: seq.chains) {
        int numOfResidues;
        file >> numOfResidues;

        chain.codes = vector<AminoAcid>(numOfResidues);
        string codesText;
        file >> codesText;
        for (int i = 0; i < numOfResidues; ++i) {
            chain.codes[i] = (AminoAcid)codesText[i];
        }

        int numOfCMaps;
        file >> numOfCMaps;

        chain.contactMaps = vector<string>(numOfCMaps);
        for (int i = 0; i < numOfCMaps; ++i) {
            filesystem::path spPath;
            file >> spPath;
            spPath = path.parent_path() / spPath;

            auto iter = seq.contactMaps.find(spPath);
            if (iter == seq.contactMaps.end()) {
                auto spFile = ifstream(spPath);
                auto sp = cmapParser.read(spFile);

                auto newNode = make_pair(spPath, move(sp));
                iter = seq.contactMaps.emplace(move(newNode)).first;
            }

            chain.contactMaps[i] = iter->first;
        }
    }

    return seq;
}
