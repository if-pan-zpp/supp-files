#pragma once
#include "../utils/Units.hpp"
#include "../data/Primitives.hpp"

namespace mdk {
    /**
     * Standard Lennard-Jones potential.
     */
    class LennardJones {
    public:
        double r_min = 5.0*angstrom, depth = 1.0*eps;

        LennardJones() = default;
        LennardJones(double r_min, double depth):
            r_min(r_min), depth(depth) {};

        /**
         * The cutoff distance for the L-J potential. I'm not sure what is the
         * "industry standard" cutoff but at 2.5 sigma the value of the
         * potential is -0.016 epsilon, which seems fine as a cutoff. For
         * comparison, at 2 sigma the value is -0.06225 epsilon which seems too
         * high to cut off.
         * @return The cutoff distance.
         */
        inline double cutoff() const {
            return 2.5 * pow(2.0, -1.0/6.0) * r_min;
        }

        /**
         * Compute the potential energy of the force field.
         * @param norm Distance between the residues.
         * @param V Variable to add the potential to.
         * @param dV_dn Variable to add the derivative to.
         */
        inline void computeV(double norm, double& V, double& dV_dn) const {
            auto norm_inv = 1.0 / norm, s = norm_inv * r_min;
            auto s6 = s*s*s*s*s*s, s12 = s6*s6;
            V += depth * (s12 - 2.0 * s6);
            dV_dn += 12 * depth * (s6 - s12) * norm_inv;
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
