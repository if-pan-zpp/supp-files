#include "forces/Chirality.hpp"
#include "simul/Simulation.hpp"
#include "data/Chains.hpp"
using namespace mdk;

void Chirality::bind(Simulation &simulation) {
    Force::bind(simulation);
    d0_cube_inv = Scalars(state->n);
    C_nat = Scalars(state->n);

    auto& model = simulation.data<Model>();

    inRange = simulation.data<Chains>().tuples(4);
    for (int i = 0; i < model.n; ++i) {
        if (!inRange[i]) continue;

        auto i1 = i-2, i2 = i-1, i3 = i, i4 = i+1;
        auto r1 = model.residues[i1].r, r2 = model.residues[i2].r,
            r3 = model.residues[i3].r, r4 = model.residues[i4].r;

        auto r12 = r2 - r1, r23 = r3 - r2, r34 = r4 - r3;
        auto r23_x_r34 = r23.cross(r34);

        auto d0_cube_inv_i = pow(r23.norm(), -3.0);
        auto C_nat_i = r12.dot(r23_x_r34) * d0_cube_inv_i;

        d0_cube_inv[i] = d0_cube_inv_i;
        C_nat[i] = C_nat_i;
    }
}

void Chirality::asyncPart(Dynamics &dyn) {
    for (int i = 0; i < (int)inRange.size(); ++i) {
        if (!inRange[i]) continue;

        auto r1 = state->r[i-2], r2 = state->r[i-1],
            r3 = state->r[i],   r4 = state->r[i+1];
        auto r12 = r2 - r1, r23 = r3 - r2, r34 = r4 - r3;
        auto r12_x_r23 = r12.cross(r23), r12_x_r34 = r12.cross(r34),
            r23_x_r34 = r23.cross(r34);

        auto C = r12.dot(r23_x_r34) * d0_cube_inv[i];
        auto diffC = C - C_nat[i];
        dyn.V += 0.5 * e_chi * diffC * diffC;

        auto f = e_chi * diffC * d0_cube_inv[i];
        dyn.F[i-2] += f * r23_x_r34;
        dyn.F[i-1] -= f * (r12_x_r34 + r23_x_r34);
        dyn.F[i] += f * (r12_x_r23 + r12_x_r34);
        dyn.F[i+1] -= r12_x_r23;
    }
}
