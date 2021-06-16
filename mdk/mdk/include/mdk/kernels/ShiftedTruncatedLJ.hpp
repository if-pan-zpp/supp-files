#pragma once
#include "../data/Primitives.hpp"
#include "../utils/Units.hpp"

namespace mdk {
    /**
     * A shifted and truncated version of the L-J potential, used chiefly
     * by the Pauli exclusion force field.
     */
    class ShiftedTruncatedLJ {
    public:
        double r_cut = 4.0 * angstrom, depth = 1.0 * eps;

        ShiftedTruncatedLJ()  = default;
        ShiftedTruncatedLJ(double r_cut, double depth):
            r_cut(r_cut), depth(depth) {};

        inline double cutoff() const {
            return r_cut;
        }

        /**
         * Compute the potential energy of the force field.
         * @param norm Distance between the residues.
         * @param V Variable to add the potential to.
         * @param dV_dn Variable to add the derivative to.
         */
        inline void computeV(double norm, double& V, double& dV_dn) const {
            auto norm_inv = 1.0 / norm, s = norm_inv * r_cut;
            auto s6 = s*s*s*s*s*s, s12 = s6*s6;
            V += depth * (s12 - 2.0 * s6 + 1.0);
            dV_dn +=  12 * depth * (s6 - s12) * norm_inv;
        }

        /**
         * Compute and add the L-J force between two residues. The
         * templates are here in order for us to be able to pass Eigen
         * expressions to it.
         * @tparam T1 Type of an lvalue to add the force on the first residue
         * to.
         * @tparam T2 Type of an lvalue to add the force on the second residue
         * to.
         * @param unit Normalized vector between the residues.
         * @param norm Distance between the residues.
         * @param V Variable to add the potential to.
         * @param F1 Lvalue to add the force on the first residue to.
         * @param F2 Lvalue to add the force on the second residue to.
         */
        template<typename T1, typename T2>
        inline void computeF(VRef unit, double norm, double& V,
            T1 F1, T2 F2) const {

            double dV_dn = 0.0;
            computeV(norm, V, dV_dn);
            F1 -= dV_dn * unit;
            F2 += dV_dn * unit;
        }
    };
}
