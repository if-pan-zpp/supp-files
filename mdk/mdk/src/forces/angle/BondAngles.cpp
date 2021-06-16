#include "forces/angle/BondAngles.hpp"
#include "data/Chains.hpp"
using namespace mdk;
using namespace std;

void BondAngles::bind(Simulation &simulation) {
    Force::bind(simulation);
    inRange = simulation.data<Chains>().triples;
}

void BondAngles::asyncPart(Dynamics &dyn) {
    #pragma omp for nowait
    for (int i = 0; i < (int) inRange.size(); ++i) {
        if (!inRange[i]) continue;

        auto r1 = state->r[i-1], r2 = state->r[i], r3 = state->r[i+1];
        auto r12 = r2 - r1, r23 = r3 - r2;

        auto r12_x_r23 = r12.cross(r23);
        double r12_x_r23_norm = r12_x_r23.norm();
        if (r12_x_r23_norm != 0.0) {
            double r12_norm = r12.norm(), r23_norm = r23.norm();

            Vector dtheta_dr1 = r12.cross(r12_x_r23).normalized() / r12_norm;
            Vector dtheta_dr3 = r23.cross(r12_x_r23).normalized() / r23_norm;
            Vector dtheta_dr2 = -dtheta_dr1 - dtheta_dr3;

            double cos_theta = -r12.dot(r23) / r12_norm / r23_norm;
            cos_theta = max(min(cos_theta, 1.0), -1.0);
            double theta = acos(cos_theta), dV_dtheta = 0.0;

            if (natBA && natBA->isNative[i]) {
                natBA->term(i, theta, dyn.V, dV_dtheta);
            }
            else if (heurBA) {
                heurBA->term(i, theta, dyn.V, dV_dtheta);
            }

            dyn.F[i-1] -= dV_dtheta * dtheta_dr1;
            dyn.F[i] -= dV_dtheta * dtheta_dr2;
            dyn.F[i+1] -= dV_dtheta * dtheta_dr3;
        }
    }
}
