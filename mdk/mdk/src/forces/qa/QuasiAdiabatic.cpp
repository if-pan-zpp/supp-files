#include "forces/qa/QuasiAdiabatic.hpp"
#include <Eigen/Core>
#include <algorithm>
using namespace mdk;
using namespace mdk::param;

void QuasiAdiabatic::bind(Simulation &simulation) {
    NonlocalForce::bind(simulation);

    stats = &simulation.var<Stats>();

    types = &simulation.data<Types>();
    chains = &simulation.data<Chains>();

    bb_lj = LennardJones(5.0 * angstrom, 1.0 * eps);
    bs_lj = LennardJones(6.8 * angstrom, 1.0 * eps);

    auto& params = simulation.data<param::Parameters>();
    for (auto const& [acids, r]: params.pairwiseMinDist) {
        auto acid1 = (int8_t)acids.first, acid2 = (int8_t)acids.second;
        auto sslj = SidechainLJ(1.0 * eps, r);
        ss_ljs[acid1][acid2] = ss_ljs[acid2][acid1] = sslj;
    }

    n = h = Vectors(state->n);

    formationMaxDistSq = 0.0;
    formationMaxDistSq = std::max(formationMaxDistSq, bb_lj.r_min);
    formationMaxDistSq = std::max(formationMaxDistSq, bs_lj.r_min);
    for (int8_t acid1 = 0; acid1 < AminoAcid::N; ++acid1) {
        for (int8_t acid2 = 0; acid2 < AminoAcid::N; ++acid2) {
            auto ss_r_min = ss_ljs[acid1][acid2].sink_max;
            formationMaxDistSq = std::max(formationMaxDistSq, ss_r_min);
        }
    }
    formationMaxDistSq = pow(formationMaxDistSq, 2.0);

    installIntoVL();
}

void QuasiAdiabatic::asyncPart(Dynamics &dyn) {
    computeNH();

    for (auto& cont: pairs) {
        if (cont.status == QAContact::Status::REMOVED)
            continue;

        double stage;
        if (cont.status == QAContact::Status::FORMING) {
            stage = std::min((state->t - cont.t0) / formationTime, 1.0);
        }
        else {
            stage = std::max(1.0 - (state->t - cont.t0) / breakingTime, 0.0);
        }

        Vector r = state->top(state->r[cont.i2] - state->r[cont.i1]);
        auto norm = r.norm();
        auto unit = r / norm;
        double r_min;

        if (stage > 0.0) {
            if (cont.type == Stats::Type::BB) {
                bb_lj.computeF(unit, norm, dyn.V, dyn.F[cont.i1],
                    dyn.F[cont.i2]);
                r_min = bb_lj.r_min;
            }
            else if (cont.type != Stats::Type::SS) {
                bs_lj.computeF(unit, norm, dyn.V, dyn.F[cont.i1],
                    dyn.F[cont.i2]);
                r_min = bs_lj.r_min;
            }
            else {
                auto const& ss_lj = ss_ljs[(*types)[cont.i1]][(*types)[cont.i2]];
                ss_lj.computeF(unit, norm, dyn.V, dyn.F[cont.i1], dyn.F[cont.i2]);
                r_min = ss_lj.sink_max;
            }

            if (cont.status == QAContact::Status::FORMING &&
                norm > breakingTolerance * pow(2.0, -1.0/6.0) * r_min) {

                cont.status = QAContact::Status::BREAKING;
                cont.t0 = state->t;
            }
        }

        if (cont.status == QAContact::Status::BREAKING && stage == 0.0) {
            cont.status = QAContact::Status::REMOVED;
            freePairs.emplace_back((QAFreePair) {
                .i1 = cont.i1, .i2 = cont.i2,
                .status = QAFreePair::Status::FREE
            });
        }
    }
}

vl::Spec QuasiAdiabatic::spec() const {
    double maxCutoff = 0.0;
    maxCutoff = std::max(maxCutoff, bb_lj.cutoff());
    maxCutoff = std::max(maxCutoff, bs_lj.cutoff());

    for (int8_t acid1 = 0; acid1 < AminoAcid::N; ++acid1) {
        for (int8_t acid2 = 0; acid2 < AminoAcid::N; ++acid2) {
            maxCutoff = std::max(maxCutoff, ss_ljs[acid1][acid2].cutoff());
        }
    }

    return (vl::Spec) {
        .cutoffSq = pow(maxCutoff, 2.0),
        .minBondSep = 3,
    };
}

static bool operator<(QAContact const& p1, std::pair<int, int> const& p2) {
    return std::make_pair(p1.i1, p1.i2) < p2;
}

static bool operator==(QAContact const& p1, std::pair<int, int> const& p2) {
    return std::make_pair(p1.i1, p1.i2) == p2;
}

void QuasiAdiabatic::vlUpdateHook() {
    std::swap(oldPairs, pairs);
    freePairs.clear();

    auto oldPairsIter = oldPairs.begin();
    auto oldPairsEnd = oldPairs.end();

    for (auto& pair : vl->pairs) {
        if (chains->isTerminal[pair.first] || chains->isTerminal[pair.second]
            || !chains->sepByAtLeastN(pair.first, pair.second, 3)) {

            continue;
        }

        while (oldPairsIter != oldPairsEnd && *oldPairsIter < pair)
            ++oldPairsIter;

        if (oldPairsIter != oldPairsEnd
            && oldPairsIter->status != QAContact::Status::REMOVED
            && *oldPairsIter == pair) {

            pairs.emplace_back(*oldPairsIter);
        }
        else {
            freePairs.emplace_back((QAFreePair) {
                .i1 = pair.first, .i2 = pair.second,
                .status = QAFreePair::Status::FREE
            });
        }
    }
}

bool QuasiAdiabatic::geometryPhase(vl::PairInfo const& p, QADiff &diff) const {
    Vector h1 = h[p.i1], h2 = h[p.i2];
    double cos_h1_r12 = h1.dot(p.unit), cos_h2_r12 = h2.dot(p.unit),
        cos_h1_h2 = h1.dot(h2);

    if (abs(cos_h1_r12) >= hr_abs_min && abs(cos_h2_r12) >= hr_abs_min &&
        abs(cos_h1_h2) >= hh_abs_min && p.norm <= bb_lj.r_min) {

        diff.cont.type = Stats::Type::BB;
        return true;
    }

    Vector n1 = n[p.i1];
    double cos_n1_r12 = n1.dot(p.unit);
    if (cos_n1_r12 <= nr_max && abs(cos_h1_r12) >= hr_abs_min &&
        p.norm <= bs_lj.r_min * formationTolerance) {

        diff.cont.type = Stats::Type::BS;
        return true;
    }

    Vector n2 = n[p.i2];
    double cos_n2_r12 = n2.dot(p.unit);
    if (cos_n2_r12 <= nr_max && abs(cos_h1_r12) >= hr_abs_min &&
        p.norm <= bs_lj.r_min * formationTolerance) {

        diff.cont.type = Stats::Type::SB;
        return true;
    }

    if (cos_n1_r12 <= nr_max && -cos_n2_r12 <= nr_max &&
        p.norm <= ss_ljs[(*types)[p.i1]][(*types)[p.i2]].sink_max * formationTolerance) {

        diff.cont.type = Stats::Type::SS;
        return true;
    }

    return false;
}

void QuasiAdiabatic::syncPart(Dynamics &dyn) {
    qaDiffs.clear();
    for (int i = 0; i < (int)freePairs.size(); ++i) {
        auto const& p = freePairs[i];

        if (p.status == QAFreePair::Status::TAKEN)
            continue;

        vl::PairInfo pairInfo;
        pairInfo.i1 = p.i1;
        pairInfo.i2 = p.i2;

        auto r = state->top(state->r[p.i1] - state->r[p.i2]);
        auto r_normsq = r.squaredNorm();
        if (r_normsq >= formationMaxDistSq)
            continue;

        pairInfo.norm = sqrt(r_normsq);
        pairInfo.unit = r / pairInfo.norm;

        QADiff diff;
        if (!geometryPhase(pairInfo, diff))
            continue;

        stats->creationDiffs(p.i1, p.i2, diff.cont.type, diff.statDiffs);

        auto stat1 = stats->stats[p.i1] + diff.statDiffs[0];
        if (!stat1.valid()) continue;

        auto stat2 = stats->stats[p.i2] + diff.statDiffs[1];
        if (!stat2.valid()) continue;

        diff.oldIdx = i;
        diff.cont.i1 = p.i1;
        diff.cont.i2 = p.i2;
        diff.cont.t0 = state->t;
        diff.cont.status = QAContact::Status::FORMING;

        qaDiffsMutex.lock();
        qaDiffs.emplace_back(diff);
        qaDiffsMutex.unlock();
    }

    std::sort(qaDiffs.begin(), qaDiffs.end());
    for (auto const& diff: qaDiffs) {
        auto& stat1 = stats->stats[diff.cont.i1];
        auto res1 = stat1 + diff.statDiffs[0];
        if (!res1.valid()) continue;

        auto& stat2 = stats->stats[diff.cont.i2];
        auto res2 = stat2 + diff.statDiffs[1];
        if (!res2.valid()) continue;

        freePairs[diff.oldIdx].status = QAFreePair::Status::TAKEN;
        stat1 = res1;
        stat2 = res2;
        pairs.push_back(diff.cont);
    }
}

void QuasiAdiabatic::computeNH() {
    for (int i = 0; i < (int)chains->triples.size(); ++i) {
        if (!chains->triples[i]) continue;

        auto r0 = state->r[i-1], r1 = state->r[i], r2 = state->r[i+1];
        auto v0 = r1 - r0, v1 = r2 - r1;
        n[i] = (v1 - v0).normalized();
        h[i] = (v1.cross(v0)).normalized();
    }
}
