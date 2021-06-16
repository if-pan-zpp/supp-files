#include "forces/go/NativeContacts.hpp"
#include "kernels/LennardJones.hpp"
using namespace mdk;

void NativeContacts::bind(Simulation &simulation) {
    NonlocalForce::bind(simulation);

    auto& model = simulation.data<Model>();
    for (auto const& cont: model.contacts) {
        if ((ContactTypeIdx)cont.type != ContactTypeIdx::SSBOND) {
            int i1 = cont.res[0];
            int i2 = cont.res[1];
            if (i1 > i2) std::swap(i1, i2);
            allContacts.emplace_back((Contact) {
                .i1 = i1, .i2 = i2,
                .r_min = cont.dist0
            });
        }
    }

    /* Remove duplicated contacts, this probably should be done in 'coarsen' */ 
    sort(allContacts.begin(), allContacts.end(),
         [](Contact const& a, Contact const& b) -> bool {
             if (a.i1 == b.i1) return a.i2 < b.i2;
             return a.i1 < b.i1;
         });
    allContacts.resize(
        distance(allContacts.begin(),
                 unique(allContacts.begin(), allContacts.end(),
                        [](Contact const& a, Contact const& b) -> bool {
                            return a.i1 == b.i1 && a.i2 == b.i2;
                        })));

    installIntoVL();
}

vl::Spec NativeContacts::spec() const {
    double maxCutoff = 18.0 * angstrom;

    return (vl::Spec) {
        .cutoffSq = pow(maxCutoff, 2.0),
        .minBondSep = 3
    };
}

bool operator<(NativeContacts::Contact const& p1, std::pair<int, int> const& p2)  {
    return std::make_pair(p1.i1, p1.i2) < p2;
}

bool operator==(NativeContacts::Contact const& p1, std::pair<int, int> const& p2) {
    return std::make_pair(p1.i1, p1.i2) == p2;
}

void NativeContacts::vlUpdateHook() {
    curPairs.clear();
    newVL.clear();

    auto allContIter = allContacts.begin();
    auto allContEnd = allContacts.end();

    for (auto const& p: vl->pairs) {
        while (allContIter != allContEnd && *allContIter < p)
            ++allContIter;

        if (allContIter != allContEnd && *allContIter == p) {
            curPairs.emplace_back(*allContIter);
        }
        else {
            newVL.emplace_back(p);
        }
    }

    std::swap(vl->pairs, newVL);
}

void NativeContacts::asyncPart(Dynamics &dyn) {
    #pragma omp for nowait 
    for (auto const& cont: curPairs) {
        auto r12 = state->top(state->r[cont.i1] - state->r[cont.i2]);
        auto x2 = r12.squaredNorm();
        if (x2 > savedSpec.cutoffSq) continue;

        auto x = sqrt(x2);
        auto unit = r12/x;

        auto lj = LennardJones(cont.r_min, depth);
        lj.computeF(unit, x, dyn.V, dyn.F[cont.i1], dyn.F[cont.i2]);
    }
}
