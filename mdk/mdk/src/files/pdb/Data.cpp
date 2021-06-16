#include "files/pdb/Data.hpp"
#include <unordered_map>
#include <unordered_set>
#include "utils/TupleHash.hpp"
using namespace mdk::pdb::records;
using namespace mdk::pdb;
using namespace mdk;
using namespace std;

pdb::Model Data::asModel() const {
    pdb::Model model;
    unordered_set<char> terminated;

    for (auto const& record: records) {
        if (holds_alternative<Atom>(record)) {
            auto& data = get<Atom>(record);

            /* This sequence of if's etc. adds an atom, whilst also possibly
             * adding the parent residue and its parent chain.
             * */
            if (model.chains.count(data.chainID) == 0) {
                model.addChain(data.chainID);
            }
            auto& chain = model.chains[data.chainID];
            if (terminated.count(data.chainID) > 0) continue;

            if (model.residues.count(data.residueSeqNum) == 0) {
                auto& res = model.addResidue(data.residueSeqNum, &chain);
                res.type = data.residueName;
            }
            auto& res = model.residues[data.residueSeqNum];

            if (model.atoms.count(data.serialNum) == 0) {
                auto& atom = model.addAtom(data.serialNum, &res);
                atom.type = data.atomName;
                atom.r = data.pos;
            }
        }
        else if (holds_alternative<Ter>(record)) {
            auto& data = get<Ter>(record);
            auto& res = model.residues[data.residueSeqNum];
            assert(data.residueName == res.type);
            terminated.insert(data.chainId);
        }
        else if (holds_alternative<SSBond>(record)) {
            auto& data = get<SSBond>(record);

            auto& ssbond = model.addContact();
            ssbond.type = "SSBOND";
            ssbond.dist0 = data.dist0;

            /* Since in the \p pdb::Model the contacts are defined between
             * atoms, we define an SSBOND to be the contact between SG atoms of
             * CYS residues; the method doesn't necessarily check whether the
             * residues are CYS, though.
             * */
            for (int i = 0; i < 2; ++i) {
                auto& info = data.res[i];
                auto& res = model.residues[info.residueSeqNum];
                ssbond.atom[i] = res.find("SG");
            }
        }
        else if (holds_alternative<Link>(record)) {
            auto& data = get<Link>(record);

            auto& link = model.addContact();
            /* The contacts between atoms don't have a well-defined notion of
             * being, say, backbone-backbone.
             */
            link.type = "NATIVE";
            link.dist0 = data.linkLength;

            for (int i = 0; i < 2; ++i) {
                auto& info = data.res[i];
                auto &res = model.residues[info.residueSeqNum];
                link.atom[i] = res.find(info.atomName);
            }
        }
        else if (holds_alternative<Cryst1>(record)) {
            auto& data = get<Cryst1>(record);
            model.top.setCell(data.size);
        }
        else if (holds_alternative<End>(record)) {
            break;
        }
    }

    return model;
}

Data::Data(const pdb::Model &model) {
    if (model.top.cell != Vector { 0.0, 0.0, 0.0 }) {
        auto record = records.emplace_back(Cryst1());
        auto& data = get<Cryst1>(record);
        data.size = model.top.cell;
        data.angles = { 0.0, 0.0, 0.0 };
    }

    for (auto const& [chainIdx, chain]: model.chains) {
        // We need some way of finding a final ATOM record for a chain.
        pdb::Model::Atom *finalAtom = nullptr;

        for (auto const& res: chain.residues) {
            for (auto const& atom: res->atoms) {
                auto& record = records.emplace_back(Atom());
                auto& data = get<Atom>(record);
                data.atomName = atom->type;
                data.residueName = atom->res->type;
                data.residueSeqNum = atom->res->serial;
                data.serialNum = atom->serial;
                data.pos = atom->r;
                data.chainID = atom->res->chain->serial;
                data.element = atom->type[0];

                finalAtom = atom;
            }
        }

        auto& record = records.emplace_back(Ter());
        auto& data = get<Ter>(record);
        auto const& finalRes = finalAtom->res;
        data.serialNum = finalAtom->serial + 1;
        data.residueSeqNum = finalRes->serial;
        data.residueName = finalRes->type;
        data.chainId = chain.serial;
    }

    for (auto const& cont: model.contacts) {
        if (cont.type == "SSBOND") {
            auto& record = records.emplace_back(SSBond());
            auto& data = get<SSBond>(record);
            data.dist0 = cont.dist0;
            data.serialNum = cont.idx;
            for (int i = 0; i < 2; ++i) {
                data.res[i].chainId = cont.atom[i]->res->chain->serial;
                data.res[i].residueSeqNum = cont.atom[i]->res->serial;
            }
        }
        else {
            auto& record = records.emplace_back(Link());
            auto& data = get<Link>(record);

            data.linkLength = cont.dist0;
            for (int i = 0; i < 2; ++i) {
                auto& info = data.res[i];
                info.atomName = cont.atom[i]->type;
                info.residueName = cont.atom[i]->res->type[i];
                info.residueSeqNum = cont.atom[i]->res->serial;
                info.chainId = cont.atom[i]->res->chain->serial;
            }
        }
    }

    records.emplace_back(End());
}

Data Data::onlyAtoms(const pdb::Model &model) {
    Data data;

    int curSerial = 0;
    for (auto const& [chainIdx, chain]: model.chains) {
        pdb::Model::Atom *finalAtom = nullptr;

        for (auto const& res: chain.residues) {
            for (auto const& atom: res->atoms) {
                auto& record = data.records.emplace_back(Atom());
                auto& atomRec = get<Atom>(record);
                atomRec.atomName = atom->type;
                atomRec.residueName = atom->res->type;
                atomRec.residueSeqNum = atom->res->serial;
                atomRec.serialNum = curSerial++;
                atomRec.pos = atom->r;
                atomRec.chainID = atom->res->chain->serial;
                atomRec.element = atom->type[0];

                finalAtom = atom;
            }
        }

        auto& record = data.records.emplace_back(Ter());
        auto& terRec = get<Ter>(record);
        auto const& finalRes = finalAtom->res;
        terRec.serialNum = curSerial++;
        terRec.residueSeqNum = finalRes->serial;
        terRec.residueName = finalRes->type;
        terRec.chainId = chain.serial;
    }

    return data;
}

namespace mdk::pdb {
    Data& operator<<(Data& data, records::Record const& record) {
        data.records.push_back(record);
        return data;
    }

    Data& operator<<(Data& data, Data const& other) {
        data.records.insert(
            data.records.end(),
            other.records.begin(),
            other.records.end());
        return data;
    }
}
