#include "files/pdb/Model.hpp"
#include "utils/AminoAcid.hpp"
#include "utils/Math.hpp"
using namespace std;
using namespace mdk::pdb;

Model::Atom& Model::addAtom(int serial, Residue *res) {
    auto& atom = atoms[serial];
    atom.serial = serial;
    atom.res = res;
    if (res) {
        atom.idxInRes = res->atoms.size();
        res->atoms.push_back(&atom);
    }
    return atom;
}

Model::Residue& Model::addResidue(int serial, Chain *chain) {
    auto& res = residues[serial];
    res.serial = serial;
    res.chain = chain;
    if (chain) {
        res.idxInChain = chain->residues.size();
        chain->residues.push_back(&res);
    }
    return res;
}

Model::Atom *Model::Residue::find(std::string const& type) {
    auto iter = find_if(atoms.begin(), atoms.end(),
        [&type](auto const& x) -> auto {
            return x->type == type;
        });

    if (iter != atoms.end()) return *iter;
    else return nullptr;
}

Model::Chain& Model::addChain(char serial) {
    auto& chain = chains[serial];
    chain.serial = serial;
    return chain;
}

Model::Contact& Model::addContact() {
    auto &cont = contacts.emplace_back();
    cont.idx = contacts.size() - 1;
    return cont;
}

Model::Model(mdk::Model const& coarse) {
    unordered_map<int, int> resIdxMap;

    for (auto const& chain: coarse.chains) {
        auto& chainHere = addChain('A' + chain.idx);
        for (int resIdx = chain.start; resIdx < chain.end; ++resIdx) {
            auto& res = coarse.residues[resIdx];
            auto& resHere = addResidue(res.idx, &chainHere);
            resHere.type = (string)res.type;

            auto& caAtom = addAtom(atoms.size(), &resHere);
            caAtom.r = res.r;
            caAtom.type = "CA";

            resIdxMap[resIdx] = caAtom.idxInRes;
        }
    }

    for (auto const& contact: coarse.contacts) {
        auto& contactHere = addContact();
        contactHere.dist0 = contact.dist0;
        contactHere.type = (string)contact.type;

        for (int i = 0; i < 2; ++i) {
            contactHere.atom[i] = &atoms[resIdxMap[contact.res[i]]];
        }
    }
}

void Model::addContactsFromAtomOverlap() {
    /* First, we populate an array of amino acid residues in the \p Model,
     * i.e. "contact atoms".
     */
    std::vector<Atom*> contactAtoms;
    for (auto& [idx, atom]: atoms) {
        if (!atom.res) continue;
        if (!AminoAcid::isProper(atom.res->type)) continue;
        contactAtoms.emplace_back(&atom);
    }

    double alpha = pow(26.0/7.0, 1.0/6.0); /* cg.f:4906 */

    /* Then, we iterate over pairs of these contact atoms and check whether
     * they are in contact; we don't do any Verlet list-like optimizations
     * because it's the prep layer.
     */
    for (int idx1 = 0; idx1 < (int)contactAtoms.size(); ++idx1) {
        auto *atom1 = contactAtoms[idx1];
        auto const& info1 = AminoAcid(atom1->res->type).info()
            .atomInfo.at(atom1->type);
        string type1 = info1.inBackbone ? "B" : "S";

        for (int idx2 = idx1; idx2 < (int)contactAtoms.size(); ++idx2) {
            auto *atom2 = contactAtoms[idx2];
            auto const& info2 = AminoAcid(atom2->res->type).info()
                .atomInfo.at(atom2->type);
            string type2 = info2.inBackbone ? "B" : "S";

            /* If the two atoms are in the same residue or in different
             * residues which are in bond distance of less than three, they
             * cannot form contact.
             */
            if (atom1->res == atom2->res)
                continue;
            if (atom1->res->chain == atom2->res->chain and
                abs(atom1->res->idxInChain - atom2->res->idxInChain) < 3)
                continue;

            /* Last, we check whether the spheres overlap, and if they do,
             * determine the type of the contact and add it.
             */
            auto dist = (atom1->r - atom2->r).norm();
            auto overlapR = info1.radius + info2.radius;
            if (dist <= alpha * overlapR) {
                auto& cont = addContact();
                cont.dist0 = dist;
                cont.type = "NAT_" + type1 + type2;
                cont.atom[0] = atom1;
                cont.atom[1] = atom2;
            }
        }
    }
}

void Model::addContactsFromResOverlap(const param::Parameters &params) {
    /* The logic is similar to \p addContactsFromAtomOverlap, but for
     * residues only.
     */
    vector<Atom*> contactAtoms;
    for (auto& [idx, res]: residues) {
        if (!AminoAcid::isProper(res.type)) continue;
        auto *caAtom = res.find("CA");
        if (caAtom) contactAtoms.emplace_back(caAtom);
    }

    for (int idx1 = 0; idx1 < (int)contactAtoms.size(); ++idx1) {
        auto *atom1 = contactAtoms[idx1];
        auto rad1 = params.radius.at((AminoAcid)atom1->type);

        for (int idx2 = idx1; idx1 < (int)contactAtoms.size(); ++idx2) {
            auto *atom2 = contactAtoms[idx2];
            auto rad2 = params.radius.at((AminoAcid)atom2->type);

            if (atom1 >= atom2) continue;
            if (atom1->res == atom2->res) continue;
            if (abs(atom1->res->idxInChain - atom2->res->idxInChain) < 3) continue;

            auto dist = (atom1->r - atom2->r).norm();
            auto overlap = rad1 + rad2;
            if (dist <= overlap) {
                auto& cont = addContact();
                cont.dist0 = dist;
                cont.type = "NAT";
                cont.atom[0] = atom1;
                cont.atom[1] = atom2;
            }
        }
    }
}

mdk::Model Model::coarsen() {
    mdk::Model model;
    model.top = top;

    std::unordered_map<int, int> resIdxMap;
    std::unordered_map<char, int> chainIdxMap;

    for (auto const& [chainIdx, chain]: chains) {
        /* For each \p pdb::Model chain we add a corresponding \p Model chain.
         */
        auto& chainThere = model.addChain();
        chainIdxMap[chain.serial] = chainThere.idx;
        int n = chain.residues.size();

        /* For each \p pdb::Model residue we add a \p Model residue, with the
         * position derived from the CA atom.
         */
        for (int resIdx = 0; resIdx < n; ++resIdx) {
            auto& res = chain.residues[resIdx];
            auto& resThere = model.addResidue(&chainThere);
            resIdxMap[res->serial] = resThere.idx;
            resThere.type = ResType(res->type);
            resThere.r = res->find("CA")->r;
            resThere.nat_r = resThere.r;
//            resThere.mass = resThere.type.mass();
            resThere.mass = 1.0;
        }

        /* Finally, we manually construct a structured part corresponding to
         * the native structure in the \p pdb::Model file.
         */
        auto& sp = model.addSP();
        chainThere.structuredParts.push_back(sp.idx);
        sp.len = n;
        sp.off = 0;
        sp.angle = vector<double>(n);
        sp.dihedral = vector<double>(n);

        for (int resIdx = 0; resIdx < n; ++resIdx) {
            auto v4 = model.residues.at(resIdx).r;
            if (resIdx > 0) {
                auto v3 = model.residues.at(resIdx-1).r;
                if (resIdx > 1) {
                    auto v2 = model.residues.at(resIdx-2).r;
                    sp.angle[resIdx-1] = angle(v2, v3, v4);
                    if (resIdx > 2) {
                        auto v1 = model.residues.at(resIdx-3).r;
                        sp.dihedral[resIdx-1] = dihedral(v1, v2, v3, v4);
                    }
                }
            }
        }
    }

    /* For each \p pdb::Model contact we add a \p Model contact between the
     * parent residues of the atoms with the corresponding type.
     */
    for (auto const& cont: contacts) {
        if (!cont.atom[0]->res || !cont.atom[1]->res) continue;
        auto& contThere = model.addContact();
        contThere.type = ContactType(cont.type);

        auto ca1 = cont.atom[0]->res->find("CA")->r;
        auto ca2 = cont.atom[1]->res->find("CA")->r;
        contThere.dist0 = (ca2 - ca1).norm();

        for (int i = 0; i < 2; ++i) {
            contThere.res[i] = resIdxMap.at(cont.atom[i]->res->idxInChain);
        }
    }

    return model;
}
