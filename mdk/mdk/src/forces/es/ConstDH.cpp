#include "forces/es/ConstDH.hpp"
using namespace mdk;

vl::Spec mdk::ConstDH::spec() const {
    return (vl::Spec) {
        .cutoffSq = pow(screeningDist, 2.0),
        .minBondSep = 3,
    };
}

void ConstDH::asyncPart(Dynamics &dyn) {
    auto coeff = pow(echarge, 2.0) / (4.0 * M_PI * permittivity);

    for (auto const& p: pairs) {
        auto r12 = state->top(state->r[p.i1] - state->r[p.i2]);
        auto x2 = r12.squaredNorm();
        if (x2 > savedSpec.cutoffSq) continue;

        auto x = sqrt(x2);
        auto unit = r12/x;

        auto V_DH = coeff * p.q1_x_q2 * exp(-x/screeningDist)/x;
        dyn.V += V_DH;

        auto dV_dx = -V_DH * (1.0 + x/screeningDist)/x;
        dyn.F[p.i1] += dV_dx * unit;
        dyn.F[p.i2] -= dV_dx * unit;
    }
}
