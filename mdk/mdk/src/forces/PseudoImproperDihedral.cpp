#include "forces/PseudoImproperDihedral.hpp"
#include "simul/Simulation.hpp"
using namespace mdk;

bool LambdaPeak::supp(double psi) const {
    double s = alpha * (psi - psi0);
    return -M_PI <= s && s <= M_PI;
}

void LambdaPeak::eval(double psi, double &L, double &dL_dpsi) const {
    double s = alpha * (psi - psi0);
    if (cosineVersion) {
        L = 0.5 * cos(s) + 0.5;
        dL_dpsi = -0.5 * alpha * sin(s);
    }
    else {
        double t = abs(s/M_PI);
        auto x_inv = 1.0/(2.0*t*t-2.0*t+1.0);
        L = (t*t-2.0*t+1.0) * x_inv;
        dL_dpsi = (2.0*t*(t-1.0)) * x_inv*x_inv;
        dL_dpsi *= (s > 0 ? 1 : 0) / M_PI;
    }
}

template<typename LJ>
void perLambda(LambdaPeak const& lambda, LJ const& lj, double psi[2],
    double norm, double& A, double& B, double& C) {
    if (lambda.supp(psi[0]) && lambda.supp(psi[1])) {
        double L[2], dL_dpsi[2];
        for (int m = 0; m < 2; ++m)
            lambda.eval(psi[m], L[m], dL_dpsi[m]);

        double ljV = 0.0, dljV_dn = 0.0;
        lj.computeV(norm, ljV, dljV_dn);

        A += dL_dpsi[0] * L[1] * ljV;
        B += dL_dpsi[1] * L[0] * ljV;
        C += L[0] * L[1] * dljV_dn;
    }
}

void PseudoImproperDihedral::deriveAngles(vl::PairInfo const& pair,
    double *psi, Vector (*dpsi_dr)[6]) const {

    int idx[2] = { pair.i1, pair.i2 };
    for (int m = 0; m < 2; ++m) {
        for (int ix = 0; ix < 6; ++ix) {
            dpsi_dr[m][ix].setZero();
        }
    }

    for (int m = 0; m < 2; ++m) {
        int i1 = idx[m]-1, i2 = idx[m], i3 = idx[m];
        int loc1 = 3*m, loc2 = 3*m+1, loc3 = 3*m+2, loc4 = 3*(1-m)+1;

        Vector r1 = state->r[i1], r2 = state->r[i2], r3 = state->r[i3];
        Vector r24 = (m == 0 ? 1 : -1) * pair.norm * pair.unit;
        Vector r12 = r2 - r1, r23 = r3 - r2, r13 = r3 - r1, r14 = r12 + r24;

        Vector rij = -r23, rkj = -r13, rkl = -r14;
        Vector rm = -r12.cross(r23);
        double rm_norm = rm.norm();
        Vector rn = r13.cross(r14);
        double rn_norm = rn.norm();

        if (rm_norm == 0.0 || rn_norm == 0.0)
            continue;

        double rkj_norm = rkj.norm();
        Vector fi = rm * rkj_norm / rm_norm;
        Vector fl = - rn * rkj_norm / rn_norm;
        Vector df = (fi * rij.dot(rkj) - fl * rkl.dot(rkj)) / (rkj_norm * rkj_norm);
        Vector fj = -fi + df;
        Vector fk = -fl - df;

        dpsi_dr[m][loc1] = fi;
        dpsi_dr[m][loc2] = fj;
        dpsi_dr[m][loc3] = fk;
        dpsi_dr[m][loc4] = fl;

        double cos_psi = rm.dot(rn) / (rm_norm * rn_norm);
        psi[m] = acos(cos_psi);
        if (rij.dot(rn) < 0.0) psi[m] = -psi[m];
    }
}

void PseudoImproperDihedral::bind(Simulation &simulation) {
    NonlocalForce::bind(simulation);

    types = &simulation.data<Types>();
    seqs = &simulation.data<Chains>();

    bb_neg_lj.r_min = 6.2 * angstrom;
    bb_neg_lj.depth = 1.0 * eps;
    bb_pos_lj.r_min = 5.6 * angstrom;
    bb_pos_lj.depth = 1.0 * eps;

    bb_pos.alpha = 6.4; bb_pos.psi0 = 1.05 * rad;
    bb_neg.alpha = 6.0; bb_neg.psi0 = -1.44 * rad;
    ss.alpha = 1.2; ss.psi0 = -0.23 * rad;

    auto& params = simulation.data<param::Parameters>();
    for (auto const& type1: AminoAcid::aminoAcids()) {
        auto code1 = (int8_t)type1;
        for (auto const& type2: AminoAcid::aminoAcids()) {
            auto code2 = (int8_t)type2;

            auto& lj = ss_ljs[code1][code2];
            lj.depth = params.mjMatrix
                       ? params.mjMatrix->at({type1, type2})
                       : 1.0 * eps;
            lj.sink_max = params.pairwiseMinDist.at({type1, type2});
        }
    }

    installIntoVL();
}

vl::Spec PseudoImproperDihedral::spec() const {
    double maxCutoff = 0.0;
    maxCutoff = std::max(maxCutoff, bb_pos_lj.cutoff());
    maxCutoff = std::max(maxCutoff, bb_neg_lj.cutoff());

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

void PseudoImproperDihedral::asyncPart(Dynamics &dyn) {
    for (auto const& [i1, i2]: pairs) {
        auto r12 = state->top(state->r[i1] - state->r[i2]);
        auto r12_normsq = r12.squaredNorm();
        if (r12_normsq >= savedSpec.cutoffSq) continue;

        auto norm = sqrt(r12_normsq);
        auto unit = r12 / norm;
        vl::PairInfo pair;
        pair.i1 = i1;
        pair.i2 = i2;
        pair.norm = norm;
        pair.unit = unit;

        double psi[2];
        Vector dpsi_dr[2][6];

        deriveAngles(pair, psi, dpsi_dr);

        /* PID potential is described by a formula:
         *   \sum_i \lambda_i(\psi_{12}) \lambda_i(\psi_{21}) \phi(r_{12})
         * Thus the derivative wrt q is:
         *   \sum_i (d\lambda_i/d\psi) d\psi_{12}/dq \lambda_i(\psi_{21}) \phi(r_{12}) +
         *          \lambda_i (d\lambda_i/d\psi) d\psi_{21}/dq \phi(r_{12}) +
         *          \lambda_i(\psi_{12}) \lambda_i(\psi_{21}) d\phi/dq
         *   = A d\psi_{12}/dq + B d\psi_{21}/dq + C d\phi/dq
         */

        double A = 0.0, B = 0.0, C = 0.0;
        auto type1 = (int8_t)(*types)[i1], type2 = (int8_t)(*types)[i2];

        perLambda(bb_pos, bb_pos_lj,
            psi, norm, A, B, C);

        perLambda(bb_neg, bb_neg_lj,
            psi, norm, A, B, C);

        perLambda(ss, ss_ljs[type1][type2],
            psi, norm, A, B, C);

        int idx[6] = { i1-1, i1, i1 + 1, i2-1, i2, i2+1 };
        for (int i = 0; i < 6; ++i) {
            dyn.F[idx[i]] -= A * dpsi_dr[0][i];
            dyn.F[idx[i]] -= B * dpsi_dr[1][i];
        }
        dyn.F[i1] += C * unit;
        dyn.F[i2] -= C * unit;
    }
}

void PseudoImproperDihedral::vlUpdateHook() {
    pairs.clear();
    for (auto const& [i1, i2]: vl->pairs) {
        bool cond = seqs->sepByAtLeastN(i1, i2, 4) &&
            !seqs->isTerminal[i1] &&
            !seqs->isTerminal[i2];

        if (cond) {
            pairs.emplace_back(i1, i2);
        }
    }
}

