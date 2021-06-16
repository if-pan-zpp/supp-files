#include "files/seq/Sequence.hpp"
using namespace std;
using namespace mdk;
using namespace mdk::seq;

Model Sequence::asModel() const {
    Model model;
    std::unordered_map<std::string, int> cmapIdxMap;

    for (auto const& [cmapName, cmap]: contactMaps) {
        auto& sp = model.addContactMap(cmap);
        cmapIdxMap[cmapName] = sp.idx;
    }

    for (auto const& chain: chains) {
        auto& chainThere = model.addChain();

        for (auto const& code: chain.codes) {
            auto& resThere = model.addResidue(&chainThere);
            resThere.type = ResType((int8_t)code);
            resThere.r = {0.0, 0.0, 0.0};
            resThere.v = {0.0, 0.0, 0.0};
//            resThere.mass = resThere.type.mass();
            resThere.mass = 1.0;
        }

        for (auto const& cmapName: chain.contactMaps) {
            model.addCMapContacts(contactMaps.at(cmapName), chainThere);
        }
    }

    return model;
}
