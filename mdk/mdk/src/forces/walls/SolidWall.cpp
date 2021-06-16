#include "forces/walls/SolidWall.hpp"
using namespace mdk;

void SolidWall::asyncPart(Dynamics &dynamics) {
    static constexpr double r0 = 5.0 * angstrom;
    static constexpr double r0_inv = 1.0 / r0;
    static constexpr double r0_sq = r0 * r0;

    for (int i = 0; i < state->n; ++i) {
        Vector v = state->r[i] - wall.projection(state->r[i]);
        auto x2 = v.squaredNorm() / r0_sq;
        if (x2 > 2.0) continue;

        if (v.dot(wall.normal()) >= 0.0) v = -v;

        auto x = sqrt(x2);
        if (x < 0.5) x = 0.5;

        auto x9 = x*x*x*x*x*x*x*x*x;
        x9 = 1.0 / x9;
        dynamics.V += eps / x9;
        auto dV_dx = -9.0 / (x9 * x);
        dynamics.F[i] += dV_dx * r0_inv * v.normalized();
    }
}
